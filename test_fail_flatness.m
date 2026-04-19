clear; close all;

%% add paths
pathinfo = dictionary();

pathinfo("mosek") = "~/ksc/matlab-install/mosek/10.1/toolbox/r2017a";
pathinfo("msspoly") = "~/ksc/matlab-install/spotless";
pathinfo("sdpt3") = "~/ksc/matlab-install/SDPT3-4.0";

pathinfo("sparsesdprelax") = "./sos-sdp-conversion";
pathinfo("modules") = "./modules";
pathinfo("popeq") = "./pop-eq";

keys = pathinfo.keys;
for i = 1: length(keys)
    key = keys(i);
    addpath(genpath(pathinfo(key)));
end

%% set Mosek parameters
param = struct(); 
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-15;  % objective gap
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-15;  % primal feasibility
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-15;  % dual   feasibility
param.MSK_DPAR_INTPNT_CO_TOL_INFEAS  = 1e-15; % infeasibility test

%% generate At
n = 5;
kappa = 2;
[M, At_sdpt3, others] = generate_data_katsura(n, kappa, param);

At_sedumi = others.At_sedumi;
m = size(At_sedumi, 2);

% build up linear system for alternating projection
linear_sys.At = At_sdpt3;
[R, ~, P] = chol(At_sdpt3' * At_sdpt3);
linear_sys.R = R;
linear_sys.P = P;

%% extract extreme point
input_info.eps = 1e-8;
input_info.eps_break = 1e-4;
input_info.eps_redundant = 1e-7;
input_info.max_iter = 20;
input_info.mosek_param = param;
input_info.linear_sys = linear_sys;

if_pop = true;
if_ray = false;

if if_pop
    [M_ray, output_info] = extract_sol_pop(others.Xopt, ...
        others.SDP.sedumi.At, others.SDP.sedumi.b, others.SDP.sedumi.c, others.SDP.sedumi.K, ...
        others.obj, input_info);
    
    X = M_ray{1};
    [V, ~] = sorteig(X);
    v = V(:, 1);
    weight = X(1, 1);
    x = v(2: n+1) / v(1);

    disp(x);
end

if if_ray 
    [ray_cellarr, output_info] = extract_ray_restart(M, At_sedumi, input_info);
    
    x_cellarr = cell(1, 1);
    w_cellarr = cell(1, 1);
    for i = 1: length(ray_cellarr)
        ray = ray_cellarr{i};
        [V, ~] = sorteig(ray);
        v = V(:, 1);
        weight = ray(1, 1);
        x = v(2: n+1) / v(1);
    
        x_cellarr{i} = x;
        disp(x');
        disp(1 - x' * x);
        w_cellarr{i} = weight;
    
        % bar(eig(ray));
    
    end
end




%% helper functions
function [M, At, others] = generate_data_BQP(n, kappa, mosek_param)
    x = msspoly('x', n);
    problem.vars = {[x]};

    

    % f = (x1 - 0.1)^4;
    f = msspoly();
    g = [msspoly()];

    h = [];
    for k = 1: n-2
        h = [h; 1 - x(k)^2];
    end
    h = [h; x(2) + x(3) + x(4) - 1];
    h = [h; x(end-1)^2 + x(end)^2 - 1];
    h = [h; 1 - x(end)^4];

    problem.objective = f;
    problem.inequality = {g};
    problem.equality = {h};

    problem.relaxation_order = kappa;
    problem.rip_predecessor = [0];
    problem.regularization.expression = msspoly(1);
    problem.regularization.value = 1;
    [SDP, info] = sparse_sdp_relax(problem);

    prob = convert_sedumi2mosek(SDP.sedumi.At, SDP.sedumi.b, SDP.sedumi.c, SDP.sedumi.K);
    [~, res] = mosekopt('minimize info', prob, mosek_param);
    [Xopt, yopt, Sopt, obj] = recover_mosek_sol_blk(res, SDP.sdpt3.blk);

    M = Xopt{1};
    
    % randomly generate a moment matrix
    At = SDP.sdpt3.At{1};
    % At = At(:, [info.At_clique.moment; info.At_clique.equality]);
    At = At(:, info.At_clique.moment);
    % [At, ~] = qr(At, 'econ');
    
    % others
    fake_At = {At};
    fake_blk = SDP.sdpt3.blk;
    fake_C = SDP.sdpt3.C;
    fake_b = SDP.sdpt3.b(1:size(At, 2));
    [At_sedumi, ~, ~, ~] = SDPT3data_SEDUMIdata(fake_blk, fake_At, fake_C, fake_b);
    % others.At_sedumi = At_sedumi;
    others.SDP = SDP;
    others.x = x;
    [mom_mat_symbolic, ~] = moment_variable(x, kappa);
    others.mom_mat_symbolic = mom_mat_symbolic;
    others.obj = max(res.sol.itr.pobjval, res.sol.itr.dobjval);
    others.Xopt = Xopt;
end

function [M, At, others] = generate_data_QS(n, kappa, mosek_param)
    x = msspoly('x', n);
    problem.vars = {[x]};

    

    x1 = x(1); x2 = x(2);

    % f = (x1 - 0.1)^4;
    f = msspoly();
    g = [msspoly()];
    h = [1 - x' * x; x1 - x2^2];

    problem.objective = f;
    problem.inequality = {g};
    problem.equality = {h};

    problem.relaxation_order = kappa;
    problem.rip_predecessor = [0];
    problem.regularization.expression = msspoly(1);
    problem.regularization.value = 1;
    [SDP, info] = sparse_sdp_relax(problem);

    prob = convert_sedumi2mosek(SDP.sedumi.At, SDP.sedumi.b, SDP.sedumi.c, SDP.sedumi.K);
    [~, res] = mosekopt('minimize info', prob, mosek_param);
    [Xopt, yopt, Sopt, obj] = recover_mosek_sol_blk(res, SDP.sdpt3.blk);

    M = Xopt{1};
    
    % randomly generate a moment matrix
    At = SDP.sdpt3.At{1};
    % At = At(:, [info.At_clique.moment; info.At_clique.equality]);
    At = At(:, info.At_clique.moment);
    % [At, ~] = qr(At, 'econ');
    
    % others
    fake_At = {At};
    fake_blk = SDP.sdpt3.blk;
    fake_C = SDP.sdpt3.C;
    fake_b = SDP.sdpt3.b(1:size(At, 2));
    [At_sedumi, ~, ~, ~] = SDPT3data_SEDUMIdata(fake_blk, fake_At, fake_C, fake_b);
    others.At_sedumi = At_sedumi;
    others.SDP = SDP;
    others.x = x;
    [mom_mat_symbolic, ~] = moment_variable(x, kappa);
    others.mom_mat_symbolic = mom_mat_symbolic;
    others.obj = max(res.sol.itr.pobjval, res.sol.itr.dobjval);
    others.Xopt = Xopt;
end