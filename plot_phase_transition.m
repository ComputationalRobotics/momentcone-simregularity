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

rank_eps = 1e-7;

sample_num = 100;
kappa = 2;

n_list = 3:10;

% type_name = "eigenvalue";
type_name = "rank";

for n = n_list
    pe_list = 2: nchoosek(n+kappa, kappa);
    % pe_list = 2: ceil( nchoosek(n+kappa, kappa) * 1.5 );

    for pe = pe_list
        filepath = sprintf("phase_transition/k=%d_n=%d_pe=%d/", kappa, n, pe);
        data_filepath = "./data/" + filepath;
        fig_filepath = "./figs/" + filepath; 
        if ~exist(fig_filepath, 'dir')
            mkdir(fig_filepath);
        end

        x_list = []; % sample id
        y_list = []; % ray rank
        

        for i = 1: sample_num 
            data_filename = data_filepath + string(i) + ".mat";
            load(data_filename);

            ray_cellarr = data.ray_cellarr;
            for j = 1: length(ray_cellarr)
                ray = ray_cellarr{j};
                d = eig(ray);
                % rescale the maximum eigenvalue to 1
                d = d / max(d); d = d';
                r = nnz(abs(d) > rank_eps);

                if type_name == "eigenvalue"
                    x_list_tmp = i * ones(1, length(d));
                    y_list_tmp = log10(abs(d) + 1e-16);
                    x_list = [x_list, x_list_tmp];
                    y_list = [y_list, y_list_tmp];
                end
                
                if type_name == "rank"
                    x_list = [x_list, i];
                    y_list = [y_list, r];
                end
            end

        end

        title_name = sprintf('$n = %d,\\ d = %d,\\ s = %d$', n, kappa, pe);

        if type_name == "rank"
            make_custom_scatter_rank(x_list, y_list, fig_filepath + "rank.png", title_name);
            % make_custom_scatter_rank_compact(x_list, y_list, fig_filepath + "rank.png", title_name);
        end
        if type_name == "eigenvalue"
            make_custom_scatter_eigenvalue(x_list, y_list, fig_filepath + "eigenvalue.png", title_name);
        end
    end
    
end


%% helper functions
function make_custom_scatter_rank(x_list, y_list, fig_name, title_name)
    options.rows = 500;
    options.cols = 400;
    options.markersize = 6;
    options.linewidth = 1.2;
    options.xls = 14;
    options.yls = 12;
    options.titlels = 15;

    rows = options.rows;
    cols = options.cols;
    markersize = options.markersize;
    linewidth = options.linewidth;
    xls = options.xls;
    yls = options.yls;
    titlels = options.titlels;

    figure('Visible', 'off', 'Position', [0, 0, cols, rows]);
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

    scatter(x_list, y_list, markersize, ...
        'Marker', 'o', ...
        'MarkerEdgeColor', [0.75, 0.10, 0.18], ...
        'MarkerFaceColor', [0.88, 0.24, 0.29], ...
        'LineWidth', 0.6, ...
        'MarkerFaceAlpha', 0.45, ...
        'MarkerEdgeAlpha', 0.55);

    ymin = floor(min(y_list));
    ymax = ceil(max(y_list));
    yticks(ymin:ymax);
    xlabel('sample id', 'FontSize', xls, 'Interpreter', 'latex');
    ylabel('rank', 'FontSize', yls, 'Interpreter', 'latex');
    title(title_name, 'FontSize', titlels, 'Interpreter', 'latex');

    ax = gca;
    ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls;
    ax.LineWidth = linewidth;
    ax.GridLineWidth = 0.5;
    ax.GridColor = [0, 0, 0];
    ax.GridAlpha = 0.3;
    ax.Box = 'off';

    grid on;

    print(fig_name, '-dpng', ['-r', num2str(300)]);
    close all;

end

function make_custom_scatter_rank_compact(x_list, y_list, fig_name, title_name)
    options.rows = 250;
    options.cols = 250;
    options.markersize = 6;
    options.linewidth = 1.2;
    options.xls = 15;
    options.yls = 15;
    options.titlels = 15;

    rows = options.rows;
    cols = options.cols;
    markersize = options.markersize;
    linewidth = options.linewidth;
    xls = options.xls;
    yls = options.yls;
    titlels = options.titlels;

    figure('Visible', 'off', 'Position', [0, 0, cols, rows]);
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

    scatter(x_list, y_list, markersize, ...
        'Marker', 'o', ...
        'MarkerEdgeColor', [0.75, 0.10, 0.18], ...
        'MarkerFaceColor', [0.88, 0.24, 0.29], ...
        'LineWidth', 0.6, ...
        'MarkerFaceAlpha', 0.45, ...
        'MarkerEdgeAlpha', 0.55);

    ymin = floor(min(y_list));
    ymax = ceil(max(y_list));
    yticks(ymin:ymax);
    xlabel('sample id', 'FontSize', xls, 'Interpreter', 'latex');
    ylabel('rank', 'FontSize', yls, 'Interpreter', 'latex');
    title(title_name, 'FontSize', titlels, 'Interpreter', 'latex');

    ax = gca;
    ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls;
    ax.LineWidth = linewidth;
    ax.GridLineWidth = 0.5;
    ax.GridColor = [0, 0, 0];
    ax.GridAlpha = 0.3;
    ax.Box = 'off';

    grid on;

    print(fig_name, '-dpng', ['-r', num2str(300)]);
    close all;

end


function make_custom_scatter_eigenvalue(x_list, y_list, fig_name, title_name)
    options.rows = 500;
    options.cols = 350;
    options.markersize = 14;
    options.linewidth = 1.2;
    options.xls = 15;
    options.yls = 12;
    options.titlels = 15;

    rows = options.rows;
    cols = options.cols;
    markersize = options.markersize;
    linewidth = options.linewidth;
    xls = options.xls;
    yls = options.yls;
    titlels = options.titlels;

    figure('Visible', 'off', 'Position', [0, 0, cols, rows]);
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

    scatter(x_list, y_list, markersize, ...
        'Marker', 'o', ...
        'MarkerEdgeColor', [0.07, 0.34, 0.61], ...
        'MarkerFaceColor', [0.18, 0.49, 0.76], ...
        'LineWidth', 0.6, ...
        'MarkerFaceAlpha', 0.40, ...
        'MarkerEdgeAlpha', 0.55);

    xlabel('sample id', 'FontSize', xls, 'Interpreter', 'latex');
    ylabel('$\log_{10}(d)$', 'FontSize', yls, 'Interpreter', 'latex');
    title(title_name, 'FontSize', titlels, 'Interpreter', 'latex');

    ax = gca;
    ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls;
    ax.LineWidth = linewidth;
    ax.GridLineWidth = 0.5;
    ax.GridColor = [0, 0, 0];
    ax.GridAlpha = 0.3;
    ax.Box = 'off';

    grid on;

    print(fig_name, '-dpng', ['-r', num2str(300)]);
    close all;

end
