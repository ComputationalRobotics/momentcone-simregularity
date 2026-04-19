clear; close all;

%% add paths
pathinfo = dictionary();

pathinfo("mosek") = "~/ksc/matlab-install/mosek/10.1/toolbox/r2017a";
pathinfo("msspoly") = "~/ksc/matlab-install/spotless";
pathinfo("sdpt3") = "~/ksc/matlab-install/SDPT3-4.0";

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

%% load data
load("./data/contact/Xs.mat");
X = Xs{1};
d = eig(X);
M = X / max(d);

%% generate moment constraints
n = 13;
kappa = 2;
mat_size = nchoosek(n + kappa, kappa);
[At_sdpt3, others] = generate_moment_cone(n, kappa, false);
At_sedumi = others.At_sedumi;
m = size(At_sedumi, 2);

% build up linear system for alternating projection
linear_sys.At = At_sdpt3;
[R, ~, P] = chol(At_sdpt3' * At_sdpt3);
linear_sys.R = R;
linear_sys.P = P;

% round M, let it become strictly feasible
[M, ~] = conic_alternating_projection(M, linear_sys, 500, 1e-14);

%% extract extreme rays
input_info.eps = 1e-7;
input_info.eps_break = 1e-4;
input_info.eps_redundant = 1e-7;
input_info.max_iter = 20;
input_info.mosek_param = param;
input_info.linear_sys = linear_sys;
[ray_cellarr, output_info] = extract_ray_restart(M, At_sedumi, input_info);

%% examine extreme rays
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
    
    % bar(log(abs(eig(ray)) + 1e-16)); 
    bar(eig(ray));

    % disp(1 - x' * x);
    % disp(x');
    
end


