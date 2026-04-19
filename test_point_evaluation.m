clear all; close all; clc;

%% add paths
pathinfo = dictionary();

% please change to your installation path of mosek and msspoly
pathinfo("mosek") = "~/ksc/matlab-install/mosek/10.1/toolbox/r2017a";
pathinfo("msspoly") = "~/ksc/matlab-install/spotless";

pathinfo("sparsesdprelax") = "./sos-sdp-conversion";
pathinfo("modules") = "./modules";

keys = pathinfo.keys;
for i = 1: length(keys)
    key = keys(i);
    addpath(genpath(pathinfo(key)));
end

%% set Mosek parameters
param = struct(); 
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-9;  % objective gap
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-15;  % primal feasibility
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-13;  % dual   feasibility
param.MSK_DPAR_INTPNT_CO_TOL_INFEAS  = 1e-15; % infeasibility test

%% generate moment constraints
n = 5;
kappa = 2;
if_sos_sdp_conversion = false; % manually set true/false
mat_size = nchoosek(n + kappa, kappa);
builder_fns = {
    @() load_constraint_cache(kappa, n), ...
    @() generate_moment_cone(n, kappa, false)
};
[At_sdpt3, others] = builder_fns{1 + if_sos_sdp_conversion}();

At_sedumi = others.At_sedumi;
m = size(At_sedumi, 2);

% build up linear system for alternating projection
linear_sys.At = At_sdpt3;
[R, ~, P] = chol(At_sdpt3' * At_sdpt3);
linear_sys.R = R;
linear_sys.P = P;

%% generate moment matrix 
point_evaluation_num = nchoosek(n+kappa, kappa);
M = zeros(mat_size);
xtrue_cellarr = cell(1, 1);
wtrue_cellarr = cell(1, 1);
for i = 1: point_evaluation_num
    weight = rand(1);
    x_value = 2 * rand(n, 1) - 1; % x is scaled in a hyper-cube for numerical stability
    point_evaluation = double( subs(others.mom_mat_symbolic, others.x, x_value) );
    M = M + weight * point_evaluation;

    xtrue_cellarr{i} = x_value;
    wtrue_cellarr{i} = weight;
end
M = 0.5 * (M + M');

%% extract extreme rays
tic;
input_info.eps = 1e-7;
input_info.eps_break = 1e-4;
input_info.eps_redundant = 1e-7;
input_info.max_iter = 20;
input_info.mosek_param = param;
input_info.linear_sys = linear_sys;
[ray_cellarr, output_info] = extract_ray_restart(M, At_sedumi, input_info);
toc;

%% compare to the true results
x_cellarr = cell(1, 1);
w_cellarr = cell(1, 1);
for i = 1: length(ray_cellarr)
    ray = ray_cellarr{i};
    [V, ~] = sorteig(ray);
    v = V(:, 1);
    weight = ray(1, 1);
    x = v(2: n+1) / v(1);

    x_cellarr{i} = x;
    w_cellarr{i} = weight;
end

x_cellarr = x_cellarr';
w_arr = [w_cellarr{:}];
[w_arr, perm] = sort(w_arr); 
w_arr = w_arr';
x_cellarr = x_cellarr(perm);

xtrue_cellarr = xtrue_cellarr';
wtrue_arr = [wtrue_cellarr{:}];
[wtrue_arr, perm] = sort(wtrue_arr); 
wtrue_arr = wtrue_arr';
xtrue_cellarr = xtrue_cellarr(perm);

%% save data for further usage
filename = sprintf("pe_n=%d_d=%d_num=%d.mat", n, kappa, point_evaluation_num);
data.ray_cellarr = ray_cellarr;
data.n = n;
data.kappa = kappa;
data.point_evaluation_num = point_evaluation_num;
data.x_cellarr = x_cellarr;
data.xtrue_cellarr = xtrue_cellarr;
data.w_arr = w_arr;
data.wtrue = wtrue_arr;
if ~exist("./data/debug/", 'dir')
    mkdir("./data/debug/");
end
save("./data/debug/" + filename, "data");

function [At_sdpt3, others] = load_constraint_cache(kappa, n)
    cone_filepath = sprintf("./constraint/moment_cone_k=%d_n=%d.mat", kappa, n);
    if ~exist(cone_filepath, 'file')
        error("Pre-stored constraint file not found: %s", cone_filepath);
    end
    cone_data = load(cone_filepath);
    if isfield(cone_data, "data")
        cone_data = cone_data.data;
    end
    At_sdpt3 = cone_data.At_sdpt3;
    others = cone_data.others;
end
















