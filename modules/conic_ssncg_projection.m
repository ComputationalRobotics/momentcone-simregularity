function [X_new, if_succeed] = conic_ssncg_projection(X, linear_sys, max_iter, pinf_eps)
    % Orthogonal projection of X onto the spectrahedral cone {X psd : A X = 0} by a
    % semismooth Newton-CG method (SSNCG) on the dual problem
    %   min_y phi(y) = 0.5 * || Pi_+(X + At * y) ||_F^2,  grad phi(y) = A * Pi_+(X + At * y).
    % The output Pi_+(X + At * y) is PSD by construction, and its primal infeasibility
    % || A * X_new || is exactly the gradient norm.
    % Same input/output as conic_alternating_projection: At is in SDPT3 svec format, and
    % R, P (sparse Cholesky factorization of At' * At) serve as the CG preconditioner.
    At = linear_sys.At;
    A = At';
    R = linear_sys.R;
    P = linear_sys.P;

    % SSNCG parameters
    tau1 = 1e-2; tau2 = 1e-2;   % regularization eps_reg = tau1 * min(tau2, ||grad||)
    eta_bar = 1e-2; tau = 0.5;  % CG relative tolerance min(eta_bar, ||grad||^tau)
    cg_maxiter = 200;
    mu = 1e-4; delta = 0.5;     % Armijo line search
    ls_maxiter = 30;
    stall_maxiter = 10;         % stop if pinf has not dropped by 10% for this many Newton iterations

    Z = svec_single(X);
    y = zeros(size(At, 2), 1);
    [Q, lam] = eig_desc(smat_single(Z));
    if_succeed = false;
    newton_iter = 0;
    cg_iter_total = 0;
    pinf_best = Inf;
    stall = 0;

    for k = 1: max_iter
        r = nnz(lam > 0);
        X_proj = Q(:, 1:r) * diag(lam(1:r)) * Q(:, 1:r)';
        grad = A * svec_single(X_proj);
        pinf = norm(grad);
        if pinf < pinf_eps
            if_succeed = true;
            break;
        end
        if pinf < 0.9 * pinf_best
            pinf_best = pinf;
            stall = 0;
        else
            stall = stall + 1;
            if stall >= stall_maxiter
                break;
            end
        end
        newton_iter = k;
        phi = 0.5 * sum(lam(1:r).^2);

        % Newton direction: (A * V * At + eps_reg * I) d = -grad, V in the generalized Jacobian of Pi_+
        Q1 = Q(:, 1:r);
        Q2 = Q(:, r+1:end);
        nu = lam(1:r) ./ (lam(1:r) - lam(r+1:end)');
        eps_reg = tau1 * min(tau2, pinf);
        Vfun = @(d) A * svec_single(psd_proj_jacobian(smat_single(At * d), Q1, Q2, nu)) + eps_reg * d;
        Mfun = @(v) chol_solve(v, R, P);
        [d, ~, ~, cg_iter] = pcg(Vfun, -grad, min(eta_bar, pinf^tau), cg_maxiter, Mfun);
        cg_iter_total = cg_iter_total + cg_iter;

        % Armijo line search on phi, starting from the full Newton step
        gd = grad' * d;
        if ~(gd < 0)
            break;
        end
        alpha = 1;
        if_accept = false;
        for ls = 1: ls_maxiter
            [Q_try, lam_try] = eig_desc(smat_single(Z + At * (y + alpha * d)));
            phi_try = 0.5 * sum(max(lam_try, 0).^2);
            % the slack 10 * eps * phi absorbs rounding errors in phi near convergence
            if phi_try <= phi + mu * alpha * gd + 10 * eps * phi
                if_accept = true;
                break;
            end
            alpha = delta * alpha;
        end
        if ~if_accept
            break;
        end
        y = y + alpha * d;
        Q = Q_try;
        lam = lam_try;
    end

    r = nnz(lam > 0);
    X_new = Q(:, 1:r) * diag(lam(1:r)) * Q(:, 1:r)';
    X_new = 0.5 * (X_new + X_new');

    if if_succeed
        fprintf("SSN-CG projection finishes! (Newton iters: %d, CG iters: %d) \n", newton_iter, cg_iter_total);
    else
        fprintf("SSN-CG projection fails! (Newton iters: %d, CG iters: %d, pinf: %3.2e) \n", ...
            newton_iter, cg_iter_total, norm(A * svec_single(X_new)));
    end
end

function VH = psd_proj_jacobian(H, Q1, Q2, nu)
    % V[H] = Q * (Omega .* (Q' * H * Q)) * Q', with Q = [Q1, Q2] and Omega = [1, nu; nu', 0],
    % computed in O(N^2 * min(r, N - r)) flops instead of O(N^3)
    if size(Q1, 2) <= size(Q2, 2)
        T = Q1' * H;
        G = 0.5 * (T * Q1) * Q1' + (nu .* (T * Q2)) * Q2';
        VH = Q1 * G;
        VH = VH + VH';
    else
        % complementary form: V[H] = H - Q * ((1 - Omega) .* (Q' * H * Q)) * Q'
        T = Q2' * H;
        G = 0.5 * (T * Q2) * Q2' + ((1 - nu)' .* (T * Q1)) * Q1';
        VH = Q2 * G;
        VH = H - (VH + VH');
    end
end

function [Q, lam] = eig_desc(W)
    [Q, D] = eig(0.5 * (W + W'));
    [lam, idx] = sort(diag(D), 'descend');
    Q = Q(:, idx);
end

function y = chol_solve(rhsy, R, P)
    rhsy = P' * rhsy;
    tmp = R' \ rhsy;
    tmp = R \ tmp;
    y = P * tmp;
end
