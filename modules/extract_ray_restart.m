function [ray_cellarr, output_info] = extract_ray_restart(M, At_sedumi, input_info)
    ray_cellarr = {};
    output_info_cellarr = cell(1, 1);
    full_rank_M = nnz( eig(M) > input_info.eps);
    linear_sys = input_info.linear_sys;
    
    iter = 1;
    input_info.eps_origin = input_info.eps;
    while true
        [ray_cellarr_single, output_info] = extract_ray_facial_reduction(M, At_sedumi, input_info);
        output_info_cellarr{iter} = output_info;
        
        if isempty(ray_cellarr)
            ray_cellarr = ray_cellarr_single;
        else
            if isempty(ray_cellarr_single)
                % Mosek fails to solve the first iteration: make eps smaller!
                input_info.eps = input_info.eps / 10;
            else
                ray_cellarr = [ray_cellarr, ray_cellarr_single];
                % reset eps
                input_info.eps = input_info.eps_origin;
            end
        end
    
        if output_info.if_success
            fprintf("Restarted extraction succeeds! \n");
            break;
        end
        if length(ray_cellarr) >= full_rank_M
            fprintf("Restarted extraction fails! \n");
            break;
        end
    
        M = output_info.M_res;
        [M, ~] = conic_alternating_projection(M, linear_sys, 500, 1e-14);
        iter = iter + 1;
    end
end



