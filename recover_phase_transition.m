clear; close all;

%% add paths
pathinfo = dictionary();

% pathinfo("mosek") = "~/ksc/matlab-install/mosek/10.1/toolbox/r2017a";
% pathinfo("msspoly") = "~/ksc/matlab-install/spotless";
% pathinfo("sdpt3") = "~/ksc/matlab-install/SDPT3-4.0";

pathinfo("mosek") = "~/matlab-install/mosek/10.1/toolbox/r2017a";
pathinfo("msspoly") = "~/matlab-install/spotless";
pathinfo("sdpt3") = "~/matlab-install/SDPT3-4.0";
pathinfo("utils") = "~/matlab-install/lab-code/utils";

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

for n = n_list
    pe_list = 2: nchoosek(n+kappa, kappa);
    % pe_list = 2: ceil( nchoosek(n+kappa, kappa) * 1.5 );

    fig_filepath = "./recover/" + sprintf("k=%d_n=%d/", kappa, n); 
    if ~exist(fig_filepath, 'dir')
        mkdir(fig_filepath);
    end

    w_error_list = [];
    x_error_list = [];

    for pe = pe_list
        filepath = sprintf("phase_transition/k=%d_n=%d_pe=%d/", kappa, n, pe);
        data_filepath = "./data/" + filepath;
        
        w_error = 0.0;
        x_error = 0.0;

        for i = 1: sample_num 
            data_filename = data_filepath + string(i) + ".mat";
            load(data_filename);

            w_arr = data.w_arr;
            wtrue_arr = data.wtrue_arr;

            if length(w_arr) ~= length(wtrue_arr)
                w_error = w_error + 1.0;
                x_error = x_error + 1.0;
            else
                x_cellarr = data.x_cellarr;
                xtrue_cellarr = data.xtrue_cellarr;

                % update w error and x error
                w_error_tmp = 0.0;
                x_error_tmp = 0.0;
                for j = 1: length(w_arr)
                    w = w_arr(j);
                    wtrue = wtrue_arr(j);
                    x = x_cellarr{j};
                    xtrue = xtrue_cellarr{j};

                    w_error_tmp = w_error_tmp + abs(w - wtrue) / (1 + abs(wtrue));
                    x_error_tmp = x_error_tmp + norm(x - xtrue) / (1 + norm(xtrue));
                end
                w_error_tmp = w_error_tmp / length(w_arr);
                x_error_tmp = x_error_tmp / length(w_arr);

                w_error = w_error + w_error_tmp;
                x_error = x_error + x_error_tmp;
            end
        end

        w_error = w_error / sample_num;
        x_error = x_error / sample_num;

        w_error_list = [w_error_list, w_error];
        x_error_list = [x_error_list, x_error];
    end

    title_name = sprintf("kappa: %d, n: %d, mat size: %d", kappa, n, nchoosek(n+kappa, kappa)); 

    fig_name = fig_filepath + "beautiful_error.png";
    make_beautiful_error_plot(pe_list, w_error_list, x_error_list, fig_name, title_name);

    % fig_name = fig_filepath + "error.png";
    % make_custom_error_plot(pe_list, w_error_list, x_error_list, fig_name, title_name);
end


%% helper functions
function make_beautiful_error_plot(pe_list, w_error_list, x_error_list, fig_name, ~)
    options.rows = 250;
    options.cols = 300;
    options.linewidth = 1.5;
    options.xls = 15;
    options.yls = 15;

    rows = options.rows;
    cols = options.cols;
    linewidth = options.linewidth;
    xls = options.xls;
    yls = options.yls;

    figure('Position', [0, 0, cols, rows]);
    set(groot, 'defaultAxesTickLabelInterpreter','latex');

    subplot(2,1,1);
    plot(pe_list, w_error_list, 'LineWidth', linewidth, 'Color', 'b');
    grid on;
    ylabel('$e_w$', 'FontSize', yls, 'Interpreter', 'latex');
    axis([-inf inf -inf inf]);
    ax = gca;
    ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls;
    ax.GridLineWidth = 0.5;
    ax.GridColor = [0, 0, 0];
    ax.GridAlpha = 0.3;

    subplot(2,1,2);
    plot(pe_list, x_error_list, 'LineWidth', linewidth, 'Color', 'r');
    grid on;
    xlabel('$s$', 'FontSize', xls, 'Interpreter', 'latex');
    ylabel('$e_z$', 'FontSize', yls, 'Interpreter', 'latex');
    axis([-inf inf -inf inf]);
    ax = gca;
    ax.XAxis.FontSize = xls;
    ax.YAxis.FontSize = yls;
    ax.GridLineWidth = 0.5;
    ax.GridColor = [0, 0, 0];
    ax.GridAlpha = 0.3;

    print(fig_name, '-dpng', ['-r', num2str(300)]);
    close all;
end

function make_custom_error_plot(pe_list, w_error_list, x_error_list, fig_name, title_name)

    % Create scatter plot
    figure('Visible', 'off');  % Use 'off' to suppress GUI if running in scripts
    
    % --- first subplot -------------------------------------------------------
    subplot(2,1,1);                     % 2 rows, 1 column, top tile
    plot(pe_list, w_error_list, 'LineWidth', 1.5);        % …your data…
    title('w error');        % sub-title for this axes
    xlabel('point num');
    ylabel('err');
    grid on;

    % --- second subplot ------------------------------------------------------
    subplot(2,1,2);                      % bottom tile
    plot(pe_list, x_error_list, 'LineWidth', 1.5);
    title('x error');
    xlabel('point num');
    ylabel('err');
    grid on;

    % --- overall title -------------------------------------------------------
    sgtitle(title_name);     % appears centred above both subplot

    % Save as PNG
    saveas(gcf, fig_name);
    close all;  % Close figure to free memory

end