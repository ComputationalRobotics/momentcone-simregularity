function [X_new, if_succeed] = conic_alternating_projection(X, linear_sys, max_iter, pinf_eps)
    % At is in SDPT3 svec format
    At = linear_sys.At;
    A = At';
    R = linear_sys.R;
    P = linear_sys.P;
    
    X_svec = svec_single(X);

    if_succeed = false;
    verbose = 0;

    for k = 1: max_iter
        % project to AX = 0: X = X0 - At * (A*At)^{-1} * A * X0
        y = chol_solve(A*X_svec, R, P);
        X_svec = X_svec - At * y;

        % project to PSD cone
        X_svec = proj_single(X_svec);

        pinf = norm(A * X_svec);
        if verbose == 1
            fprintf("iter: %d, pinf: %3.2e \n", k, pinf);
        end
        if pinf < pinf_eps
            if_succeed = true;
            fprintf("Alternating projection finishes! \n");
            break;
        end
    end

    if ~if_succeed 
        fprintf("Alternating projection fails! \n");
    end

    X_new = smat_single(X_svec);
end

function y = chol_solve(rhsy, R, P)
    rhsy = P' * rhsy;
    tmp = R' \ rhsy;
    tmp = R \ tmp;
    y = P * tmp;
end

