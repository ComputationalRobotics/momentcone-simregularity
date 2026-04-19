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

kappa = 2;
n_list = 3: 10;

% img_type = "eigenvalue";
img_type = "rank";

for n = n_list 
    pe_list = 2: nchoosek(n+kappa, kappa);
    % pe_list = 2: ceil( nchoosek(n+kappa, kappa) * 1.5 );
    
    % generate figure addresses
    img_paths = [];
    for pe = pe_list
        filepath = sprintf("phase_transition/k=%d_n=%d_pe=%d/", kappa, n, pe);
        img_name = "./figs/" + filepath + img_type + ".png";
        img_paths = [img_paths; img_name];
    end

    gif_path = "./videos/" + sprintf("phase_transition/k=%d_n=%d/", kappa, n);
    if ~exist(gif_path, 'dir')
        mkdir(gif_path);
    end
    gif_name = gif_path + img_type + ".gif";
    makeGif(img_paths, gif_name, 0.5); 
end




function makeGif(imgPaths, outGif, delayTime)
    %MAKEGIF  Assemble a series of still images into an animated GIF.
    %
    %   makeGif(IMGPATHS) combines the images listed in the string array or
    %   cell array IMGPATHS into "animation.gif", using a 0.05-s frame delay.
    %
    %   makeGif(IMGPATHS, OUTGIF) writes the result to OUTGIF instead.
    %
    %   makeGif(IMGPATHS, OUTGIF, DELAY) sets the inter-frame delay (seconds).
    %
    %   Example
    %   -------
    %     frames = ["figs/frame_01.png", "figs/frame_02.png", "figs/frame_03.png"];
    %     makeGif(frames, 'trajectory.gif', 0.08);
    %
    %   Notes
    %   -----
    %   • Images are written in the order given.  If you need “natural” sorting
    %     (frame_2 before frame_10) call `imgPaths = natsort(imgPaths);`
    %     beforehand, e.g. from File Exchange’s NATSORTFILES.
    %   • Transparent PNGs keep their transparency if they contain a valid
    %     alpha channel.

    % -------- argument handling --------
    arguments
        imgPaths {mustBeNonempty}
        outGif  (1,1) string    = "animation.gif"
        delayTime (1,1) double  = 0.05        % seconds per frame
    end

    % allow cellstr or string input
    imgPaths = cellstr(imgPaths);

    % -------- write the GIF --------
    for k = 1:numel(imgPaths)
        % read frame
        [A, ~, alpha] = imread(imgPaths{k});

        % if RGBA (PNG with transparency), convert to indexed w/ alpha mask
        if ~isempty(alpha)
            % merge alpha into image background—uses white (255) here
            A = im2uint8(A);  % ensure uint8
            A(repmat(alpha == 0, [1 1 3])) = 255;
        end

        % to indexed image & colormap (256 colors max for GIF)
        [I, map] = rgb2ind(A, 256, 'nodither');

        % write or append
        if k == 1
            imwrite(I, map, outGif, ...
                    'gif', ...
                    'LoopCount', Inf, ...
                    'DelayTime', delayTime);
        else
            imwrite(I, map, outGif, ...
                    'gif', ...
                    'WriteMode', 'append', ...
                    'DelayTime', delayTime);
        end
    end

    fprintf('GIF saved: %s (%d frames, %.0f ms per frame)\n', ...
            outGif, numel(imgPaths), delayTime*1000);
end
