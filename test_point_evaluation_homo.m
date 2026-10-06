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

%% choose homogeneous atoms z = (x0, x)
% "tex":    Test 2 of homogeneous_moment_minimal_test.tex, one atom at infinity
% "random": random atoms in a hyper-cube, inf_num of them at infinity (x0 = 0)
example = "tex";
if example == "tex"
    n = 2;
    kappa = 2;
    ztrue = [1, 1, 0; 0, 1, 1; 0, 0, 1]; % columns: (1,0,0), (1,1,0), (0,1,1)
    wtrue = [0.5, 0.5, 1];
else
    n = 5;
    kappa = 2;
    point_evaluation_num = n + kappa;
    % atoms at infinity all lie on the hyperplane x0 = 0: too many of them make
    % the decomposition non-unique
    inf_num = 2;
    ztrue = 2 * rand(n+1, point_evaluation_num) - 1;
    ztrue(1, 1: inf_num) = 0;
    wtrue = rand(1, point_evaluation_num);
end
atom_num = size(ztrue, 2);

% represent each atom with its largest-magnitude entry equal to 1
for i = 1: atom_num
    [~, j] = max(abs(ztrue(:, i)));
    c = ztrue(j, i);
    ztrue(:, i) = ztrue(:, i) / c;
    wtrue(i) = wtrue(i) * c^(2 * kappa);
end

%% generate moment constraints
if_sos_sdp_conversion = false; % manually set true/false
mat_size = nchoosek(n + kappa, kappa);
builder_fns = {
    @() load_constraint_cache(kappa, n), ...
    @() generate_moment_cone(n, kappa, false)
};
[At_sdpt3, others] = builder_fns{1 + if_sos_sdp_conversion}();
At_sedumi = others.At_sedumi;

% build up linear system for alternating projection
linear_sys.At = At_sdpt3;
[R, ~, P] = chol(At_sdpt3' * At_sdpt3);
linear_sys.R = R;
linear_sys.P = P;

%% generate homogeneous moment matrix
% basis of the moment matrix: x0^(kappa - |alpha|) * x^alpha, |alpha| <= kappa.
% row k of E holds the exponents of (x0, x1, ..., xn) of the k-th basis element
E = homogeneous_exponents(n, kappa);
M = zeros(mat_size);
for i = 1: atom_num
    v = homogeneous_veronese(ztrue(:, i), E);
    M = M + wtrue(i) * (v * v');
end
M = 0.5 * (M + M');

% Jean's normalization y_0 = M(1, 1) = 1 (atoms at infinity do not contribute to y_0)
y0 = M(1, 1);
assert(y0 > 0, "at least one atom must have x0 ~= 0");
M = M / y0;
wtrue = wtrue / y0;

% M must satisfy the moment constraints; otherwise E is in the wrong order
pinf = norm(At_sedumi' * M(:)) / norm(M, 'fro');
assert(pinf < 1e-10, "M violates the moment constraints: %3.2e", pinf);

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

%% compare to the true results with projective decoding
rank_eps = 1e-7;
ray_num = length(ray_cellarr);
zrec_all = zeros(n+1, ray_num);
wrec_all = zeros(1, ray_num);
chart_all = zeros(1, ray_num);
ray_rank = zeros(1, ray_num);
M_rec = zeros(mat_size);
for k = 1: ray_num
    ray = ray_cellarr{k};
    [zrec_all(:, k), wrec_all(k), chart_all(k)] = decode_projective(ray, E);
    lam = eig(ray);
    ray_rank(k) = nnz(lam > rank_eps * max(lam));
    M_rec = M_rec + ray;
end

% match recovered and planted atoms: minimum-cost assignment on e_proj
cost = zeros(ray_num, atom_num);
for k = 1: ray_num
    for i = 1: atom_num
        cost(k, i) = proj_err(zrec_all(:, k), ztrue(:, i));
    end
end
pairs = matchpairs(cost, 10);

zrec = nan(n+1, atom_num);
wrec = nan(1, atom_num);
chart = nan(1, atom_num);
e_proj = nan(1, atom_num);
e_w = nan(1, atom_num);
e_proj_affine = nan(1, atom_num);
for p = 1: size(pairs, 1)
    k = pairs(p, 1);
    i = pairs(p, 2);
    zrec(:, i) = zrec_all(:, k);
    wrec(i) = wrec_all(k);
    chart(i) = chart_all(k);
    e_proj(i) = cost(k, i);
    e_w(i) = abs(wrec(i) - wtrue(i)) / wtrue(i);
    % affine decoding of test_point_evaluation.m, for comparison
    [V, ~] = sorteig(ray_cellarr{k});
    x = V(2: n+1, 1) / V(1, 1);
    e_proj_affine(i) = proj_err([1; x], ztrue(:, i));
end

fprintf("rank(M) = %d, #rays = %d, ray ranks = [%s], reconstruction error = %3.2e \n", ...
    nnz(eig(M) > input_info.eps), ray_num, num2str(ray_rank), norm(M - M_rec, 'fro') / norm(M, 'fro'));
fprintf("atom | planted z | recovered z | chart | e_proj | e_w | e_proj (affine decoding) \n");
for i = 1: atom_num
    fprintf("%d | [%s ] | [%s ] | x%d | %3.2e | %3.2e | %3.2e \n", i, ...
        num2str(ztrue(:, i)', "%8.4f"), num2str(zrec(:, i)', "%8.4f"), chart(i) - 1, ...
        e_proj(i), e_w(i), e_proj_affine(i));
end

%% save data for further usage
filename = sprintf("pe_homo_%s_n=%d_d=%d_num=%d.mat", example, n, kappa, atom_num);
data.ray_cellarr = ray_cellarr;
data.n = n;
data.kappa = kappa;
data.ztrue = ztrue;
data.wtrue = wtrue;
data.zrec = zrec;
data.wrec = wrec;
data.e_proj = e_proj;
data.e_w = e_w;
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

function E = homogeneous_exponents(n, kappa)
    % sos-sdp-conversion stores each basis element as a sorted list of kappa
    % indices in {0, ..., n}, where index 0 is x0; this enumerates them in the
    % same order as next_choice in sos-sdp-conversion/include/indices.hpp
    mat_size = nchoosek(n + kappa, kappa);
    E = zeros(mat_size, n + 1);
    seq = zeros(1, kappa);
    for k = 1: mat_size
        E(k, :) = accumarray(seq' + 1, 1, [n + 1, 1])';
        idx = find(seq ~= n, 1, 'last');
        if ~isempty(idx)
            seq(idx: end) = seq(idx) + 1;
        end
    end
end

function v = homogeneous_veronese(z, E)
    % v(k) = z0^E(k, 1) * z1^E(k, 2) * ... * zn^E(k, n+1)
    v = prod(z(:)' .^ E, 2);
end

function [z, w, chart] = decode_projective(ray, E)
    % ray = w * v(z) * v(z)': take the largest pure-power entry h(kappa * e_j) of the
    % top eigenvector h, then z_l / z_j = h((kappa-1) * e_j + e_l) / h(kappa * e_j).
    % this never divides by the x0^kappa entry, so atoms at infinity are fine
    kappa = sum(E(1, :));
    var_num = size(E, 2);
    [V, ~] = sorteig(ray);
    h = V(:, 1);
    pure = zeros(1, var_num);
    for j = 1: var_num
        pure(j) = find(E(:, j) == kappa);
    end
    [~, chart] = max(abs(h(pure)));
    mixed = zeros(1, var_num);
    for l = 1: var_num
        e = zeros(1, var_num);
        e(chart) = kappa - 1;
        e(l) = e(l) + 1;
        [~, mixed(l)] = ismember(e, E, 'rows');
    end
    z = h(mixed) / h(pure(chart));
    % same representative as the planted atoms: largest-magnitude entry equal to 1
    [~, j] = max(abs(z));
    z = z / z(j);
    v = homogeneous_veronese(z, E);
    w = (v' * ray * v) / (v' * v)^2;
end

function e = proj_err(zhat, z)
    % min over lambda of ||zhat - lambda * z|| / ||zhat||
    if any(~isfinite(zhat))
        e = 1;
        return;
    end
    lambda = (z' * zhat) / (z' * z);
    e = norm(zhat - lambda * z) / norm(zhat);
end
