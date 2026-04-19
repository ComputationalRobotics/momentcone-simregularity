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

%% load data
n = 10;
kappa = 2;
point_evaluation_num = 66;
filename = sprintf("pe_n=%d_d=%d_num=%d.mat", n, kappa, point_evaluation_num);
moment_matrix_info = load("./data/debug/" + filename);
moment_matrix_info = moment_matrix_info.data;
ray_cellarr = moment_matrix_info.ray_cellarr;

%% generate moment constraints
mat_size = nchoosek(n + kappa, kappa);
if_sos_sdp_conversion = false; % manually set true/false
builder_fns = {
    @() load_constraint_cache(kappa, n), ...
    @() generate_moment_cone(n, kappa, true)
};
[~, others] = builder_fns{1 + if_sos_sdp_conversion}();
Bt_sdpt3 = others.Bt_sdpt3;
Bt_sedumi = others.Bt_sedumi;

%% Important: to make set compact: y_1 = 1!: remove the first constraint
Bt_sedumi = Bt_sedumi(:, 2:end);

%% get u1 ... ur spanning kernel M
for ray_num = 1: length(ray_cellarr)

    M = ray_cellarr{ray_num};
    bar(log10(abs(eig(M) + 1e-16)));

    eps = 1e-6;
    [V, D] = eig(M);
    d = diag(D);
    r = nnz(abs(d) < eps);
    U = V(:, 1:r);
    
    % generate N 
    nd = nchoosek(n+kappa, kappa);
    n2d = nchoosek(n+2*kappa, 2*kappa);
    assert(mat_size == nd);
    assert(size(Bt_sedumi, 2) == n2d - 1);
    N = zeros(r * nd, n2d-1);
    for k = 1: n2d-1 
        Bk = Bt_sedumi(:, k);
        Bk = reshape(Bk, nd, nd);
        mat = Bk * U;
        N(:, k) = mat(:);
    end
    s = svd(N);
    bar(log10(s + 1e-16));
    
    fprintf("ray %d: rank %d, sig_min: %3.2e \n", ray_num, nd - r, min(s));

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






