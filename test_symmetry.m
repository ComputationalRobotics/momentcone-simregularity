clear; close all;

%% add paths
pathinfo = dictionary();

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
mat_size = nchoosek(n + kappa, kappa);
[At_sdpt3, others] = load_constraint_cache(kappa, n);
At_sedumi = others.At_sedumi;
m = size(At_sedumi, 2);

% build up linear system for alternating projection
linear_sys.At = At_sdpt3;
[R, ~, P] = chol(At_sdpt3' * At_sdpt3);
linear_sys.R = R;
linear_sys.P = P;

%% generate M
SDP = BQP(n, kappa);
prob   = convert_sedumi2mosek(SDP.sedumi.At, SDP.sedumi.b, SDP.sedumi.c, SDP.sedumi.K);
[~, res] = mosekopt('minimize info', prob, param);
[Xopt, ~, ~, ~] = recover_mosek_sol_blk(res, SDP.sdpt3.blk);
M = Xopt{1};

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

%% check the rank and value
x_cellarr = [];
w_cellarr = [];
for i = 1: length(ray_cellarr)
    ray = ray_cellarr{i};
    [V, ~] = sorteig(ray);
    v = V(:, 1);
    weight = ray(1, 1);
    x = v(2: n+1) / v(1);

    x_cellarr = [x_cellarr, x];
    w_cellarr = [w_cellarr, weight];
end







%% helper functions
function SDP = BQP(n, kappa)
    x = msspoly('x', n);

    problem.vars = {x};
    problem.objective = msspoly(0);
    problem.inequality = {msspoly()};
    problem.equality = {x.^2 - 1};
    problem.relaxation_order = kappa;
    problem.rip_predecessor = 0;
    problem.regularization.expression = msspoly(1);
    problem.regularization.value = 1;

    [SDP, ~] = sparse_sdp_relax(problem);
    % pool = generate_pool_fast(SDP, 1e-10, 5e-4);
    % SDP = remove_redundant(SDP, pool);
    % fprintf("\n At's smallest singular value: %3.2e \n", svds(SDP.sedumi.At, 1, 'smallest'));
end

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