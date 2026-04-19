function [M_ray, output_info] = extract_sol_pop(M, At, b, c, K, obj, input_info)
    % M is of the standard matrix cell array format
    % At, b, K are of the sedumi format directly from POP.
    % Change K.s as [K.s, 1], where the last 1*1 PSD matrix is a relaxation from.
    % Change At as [ At, vec(c); 0, 1].
    % Change b as [b; obj].
    % Change M as [M, {obj - <C, X>}].
    % <C, X> <= p_obj <--> <C, X> + t = obj, t >= 0
    % Therefore, when generating a random linear functional, we need to make sure its last element is 0. 
    % Since now we just need to randomly hit one extreme point, we only need to extract once.
    
    %% change sedumi data
    % change M
    M_vec = cell2vec(M);
    obj_fake = c' * M_vec;
    if obj_fake > obj 
        warning("<C, X> > obj! \n");
        fprintf("<C, X>: %5.4e, obj: %5.4e \n", obj_fake, obj);
        obj = obj_fake + 1e-10;
        M = [M, {obj - obj_fake}];
    else
        M = [M, {obj - obj_fake}];
    end
    K.s = [K.s, 1];
    At = [
        At, c(:);
        zeros(1, size(At, 2)), 1;
    ];
    b = [b; obj];

    %% hyper parameters
    eps = input_info.eps;
    eps_break = input_info.eps_break;
    eps_redundant = input_info.eps_redundant;
    mosek_param = input_info.mosek_param;
    max_iter = input_info.max_iter;
    linear_sys = input_info.linear_sys;
    output_info.if_succeed = false;

    Q = cell(size(M));
    rank_list = zeros(size(M));
    
    for iter = 1: max_iter
        %% generate new Q from M's column spaces
        K_comp.s = K.s;
        for i = 1: length(M)
            Ms = M{i};
            if size(Ms, 1) > 1
                [V, D] = sorteig(Ms);
                d = diag(D);
                rank_Ms = nnz(d > eps);
                Qs = V(:, 1: rank_Ms);
                rank_list(i) = rank_Ms;
                K_comp.s(i) = rank_Ms;
                Q{i} = Qs;
            else
                rank_list(i) = 1;
                K_comp.s(i) = 1;
                Q{i} = 1;
            end
        end

        %% compress At_sedumi
        n_comp = 0;
        for i = 1: length(rank_list)
            n_comp = n_comp + rank_list(i)^2;
        end
        At_comp = zeros(n_comp, size(At, 2));
        for j = 1: size(At, 2)
            Aj = vec2cell(At(:, j), K);
            Aj_comp = cell(size(Q));
            for i = 1: length(M)
                Qs = Q{i};
                Ajs = Aj{i};
                Aj_comp{i} = Qs' * Ajs * Qs;
            end
            At_comp(:, j) = cell2vec(Aj_comp);
        end
        % remove redundant constraints
        [At_comp, ids] = remove_dependent_cols(At_comp, eps_redundant);
        b_comp = b(ids);

        %% solve Mosek to hit an extreme point
        % add random linear objective
        Face.c = cell(size(M));
        for i = 1: length(M)-1
            mat = randn(K_comp.s(i));
            mat = 0.5 * (mat + mat');
            Face.c{i} = mat;
        end
        Face.c{length(M)} = 1;
        Face.c = cell2vec(Face.c);
        Face.K.s = K_comp.s;
        Face.blk = cell(length(K_comp.s), 2);
        for i = 1: length(K.s)
            Face.blk{i, 1} = 's';
            Face.blk{i, 2} = K_comp.s(i);
        end
        Face.At = At_comp;
        Face.b = b_comp;
        % solve
        prob   = convert_sedumi2mosek(Face.At, Face.b, Face.c, Face.K);
        [~, res] = mosekopt('minimize info', prob, mosek_param);
        [M_ray_comp, ~, ~, ~] = recover_mosek_sol_blk(res, Face.blk);

        %% recover the solutions from compressed data
        M_ray = cell(size(M_ray_comp));
        for i = 1: length(M_ray)
            Qs = Q{i};
            M_rays_comp = M_ray_comp{i};
            M_ray{i} = Qs * M_rays_comp * Qs';
        end
        % remove the auxilary variable t
        M_ray = M_ray(1: end-1);

        %% check exist condition
        if abs(res.sol.itr.pobjval - res.sol.itr.dobjval) > eps_break
            fprintf("Extreme point generation fails since Mosek fails! \n");
            eps = eps / 10;
        else 
            fprintf("Extreme point generation succeed! \n");
            output_info.if_succeed = true;
            break;
        end
    end


end


function M_vec = cell2vec(M)
    M_vec = [];
    for i = 1: length(M)
        Ms = M{i};
        M_vec = [M_vec; Ms(:)];
    end
end

function M = vec2cell(M_vec, K)
    M = cell(size(K.s));
    idx = 1;
    mat_sizes = K.s;
    for i = 1: length(K.s)
        Ms_vec = M_vec(idx: idx + mat_sizes(i)^2 - 1, 1);
        M{i} = reshape(Ms_vec, mat_sizes(i), mat_sizes(i));
        idx = idx + mat_sizes(i)^2;
    end
end







