clearvars;
close all;

%% add paths
pathinfo = dictionary();

% Please change these to your local installation paths.
pathinfo("mosek") = "~/ksc/matlab-install/mosek/10.1/toolbox/r2017a";
pathinfo("msspoly") = "~/ksc/matlab-install/spotless";

pathinfo("sparsesdprelax") = "./sos-sdp-conversion";
pathinfo("modules") = "./modules";

keys = pathinfo.keys;
for i = 1:length(keys)
    key = keys(i);
    addpath(genpath(pathinfo(key)));
end

%% batch generation setup
kappa = 4;
n_list = 2:3;
outdir = "./constraint";
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

fprintf("Generating constraint cache for kappa=%d, n=%d:%d\n", ...
    kappa, n_list(1), n_list(end));

%% generate and save
for n = n_list
    fprintf("[start] n=%d, kappa=%d\n", n, kappa);
    tic;

    [At_sdpt3, others] = generate_moment_cone(n, kappa, true);

    % Keep only fields needed by test_point_evaluation and test_verify_extremity.
    others_keep = struct();
    others_keep.At_sedumi = others.At_sedumi;
    others_keep.Bt_sdpt3 = others.Bt_sdpt3;
    others_keep.Bt_sedumi = others.Bt_sedumi;
    others_keep.x = others.x;
    others_keep.mom_mat_symbolic = others.mom_mat_symbolic;

    data = struct();
    data.n = n;
    data.kappa = kappa;
    data.At_sdpt3 = At_sdpt3;
    data.others = others_keep;

    filename = sprintf("moment_cone_k=%d_n=%d.mat", kappa, n);
    save(fullfile(outdir, filename), "data", "-v7.3");

    fprintf("[saved] %s (%.2fs)\n", fullfile(outdir, filename), toc);
end

fprintf("All done. Files are saved in %s\n", outdir);
