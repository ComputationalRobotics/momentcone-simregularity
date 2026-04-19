function [ray_cellarr, output_info] = extract_ray_facial_reduction(M, At_sedumi, input_info)

    eps = input_info.eps;
    eps_break = input_info.eps_break;
    eps_redundant = input_info.eps_redundant;
    mosek_param = input_info.mosek_param;
    linear_sys = input_info.linear_sys;
    mat_size = size(M, 1);
    mat_size_origin = mat_size;
    Q_full = eye(mat_size);
    ray_cellarr = {};
    ids_cellarr = cell(1, 1);
    Q_cellarr = cell(1, 1);
    ray_idx = 1;
    if_success = false;
    m = size(At_sedumi, 2);
    At_sedumi_origin = At_sedumi;

    ids_global = 1: m;

    while true
        % get new rank
        [V, D] = sorteig(M);
        d = diag(D);
        rank_M = nnz(d > eps);
        if rank_M == 1
            % return M itself as the extreme ray
            ray = Q_full * M * Q_full';
            ray_cellarr{ray_idx} = 0.5 * (ray + ray');
            fprintf("Extreme ray decomposition finishes since M is already rank-1! \n");
            if_success = true;
            break;
        end
        if ray_idx == 1
            max_iter = min(input_info.max_iter, rank_M);
        end
        Q = V(:, 1: rank_M);
        Q_cellarr{ray_idx} = Q;
        Q_full = Q_full * Q;
        
        % compress At_sedumi
        At_sedumi_selected = At_sedumi_origin(:, ids_global);
        At_sedumi_compress = zeros(rank_M^2, m);
        for i = 1: m 
            A_tmp = At_sedumi_selected(:, i);
            A_tmp = reshape(A_tmp, mat_size_origin, mat_size_origin);
            A_tmp_compress = Q_full' * A_tmp * Q_full;
            A_tmp_compress = 0.5 * (A_tmp_compress + A_tmp_compress');
            At_sedumi_compress(:, i) = A_tmp_compress(:);
        end

        % remove dependent constraints
        [At_sedumi_compress, ids] = remove_dependent_cols(At_sedumi_compress, eps_redundant);
        ids_cellarr{ray_idx} = ids;
        ids_global_tmp = zeros(1, length(ids));
        for ii_ids = 1: length(ids)
            ids_global_tmp(ii_ids) = ids_global( ids(ii_ids) );
        end
        ids_global = ids_global_tmp;
        
        % add trace constraint
        M_compress = Q' * M * Q;
        M_compress = 0.5 * (M_compress + M_compress');
        I = eye(rank_M);
        Face.At = [At_sedumi_compress, I(:)];
        Face.b = zeros(size(Face.At, 2), 1);
        Face.b(end) = trace(M_compress);
        
        % add random linear objective
        Face.c = randn(rank_M);
        Face.c = 0.5 * (Face.c + Face.c');
        Face.c = Face.c(:);
        Face.K.s = rank_M;
        Face.blk = cell(1, 2);
        Face.blk{1, 1} = 's';
        Face.blk{1, 2} = rank_M;
        
        % solve Mosek to hit an extreme ray
        prob   = convert_sedumi2mosek(Face.At, Face.b, Face.c, Face.K);
        [~, res] = mosekopt('minimize info', prob, mosek_param);
        [Xopt, ~, ~, ~] = recover_mosek_sol_blk(res, Face.blk);
        ray_compress = Xopt{1};
        
        % hold off;
        % semilogy(abs(eig(ray_compress)) + 1e-16); hold on; grid on;

        % improve ray_compress's feasibility
        ray = Q_full * ray_compress * Q_full';
        ray = 0.5 * (ray + ray');
        err = norm(At_sedumi_origin' * ray(:));
        if ( err > eps_break ) || ... 
           ( abs(res.sol.itr.pobjval - res.sol.itr.dobjval) > eps_break ) 
            fprintf("Extreme ray decomposition fails since Mosek fails! \n");
            M_full = Q_full * M_compress * Q_full';
            M_full = 0.5 * (M_full + M_full');
            break;
        else
            input_info.eps = input_info.eps_origin;
        end

        [ray_clean, ~] = conic_alternating_projection(ray, linear_sys, 500, 1e-14);
        ray_compress = Q_full' * ray_clean * Q_full;
        ray_compress = 0.5 * (ray_compress + ray_compress');

        % semilogy(abs(eig(ray_compress)) + 1e-16); 
    
        % if res.info.MSK_DINF_INTPNT_PRIMAL_FEAS > eps
        %     % return M itself as the "unfinished" extreme ray
        %     ray = Q_full * M_compress * Q_full';
        %     ray_cellarr{ray_idx} = 0.5 * (ray + ray');
        %     fprintf("Extreme ray decomposition fails! \n");
        %     break;
        % end
    
        if norm(M_compress - ray_compress, 'fro') < eps_break
            % M_compress is already an extreme ray
            ray = Q_full * ray_compress * Q_full';
            ray_cellarr{ray_idx} = 0.5 * (ray + ray');
            fprintf("Extreme ray decomposition finishes since M is already an extreme ray! \n");
            if_success = true;
            break;
        end
        
        % extract the extreme ray for M_compress:
        % alpha = lam_max^{-1} ( M_compress^(-0.5) * ray_compress * M_compress^(-0.5) )
        [V, D] = sorteig(M_compress);
        d_sqrt = sqrt( diag(D) );
        D_sqrtinv = diag( 1 ./ d_sqrt );
        M_compress_sqrtinv = V * D_sqrtinv * V';
        mat = M_compress_sqrtinv * ray_compress * M_compress_sqrtinv;
        mat = 0.5 * (mat + mat');
        lam = eig(mat);
        lammax = lam(end);
        ray_compress = 1 / lammax * ray_compress;
        ray = Q_full * ray_compress * Q_full';
        ray_cellarr{ray_idx} = 0.5 * (ray + ray');
        M_compress_new = M_compress - ray_compress;
        
        % remove small eigenvalues
        [V, D] = sorteig(M_compress_new);
        d = abs(diag(D));
        % semilogy(d);
        d(d < eps) = 0;
        M_compress_new = V * diag(d) * V';
        M_compress_new = 0.5 * (M_compress_new + M_compress_new');
        
        % update information
        At_sedumi = At_sedumi_compress;
        M = M_compress_new;
        mat_size = rank_M;
        m = size(At_sedumi, 2);
        ray_idx = ray_idx + 1;
        
        % for debug
        M_full = Q_full * M * Q_full';
        M_full = 0.5 * (M_full + M_full');
        err = norm(At_sedumi_origin' * M_full(:));
        fprintf("primal infeasibility: %3.2e \n", err);
        % semilogy(eig(M_full));
    
        if ( err > eps_break ) || ...
           ( ray_idx > max_iter )
            fprintf("Extreme ray decomposition unfinishes due to restart! \n");
            break;
        end
    end

    output_info.if_success = if_success;
    output_info.Q_cellarr = Q_cellarr;
    output_info.ids_cellarr = ids_cellarr;
    if if_success
        output_info.M_res = zeros(size(mat_size_origin, 1));
    else
        output_info.M_res = M_full;
    end
end

