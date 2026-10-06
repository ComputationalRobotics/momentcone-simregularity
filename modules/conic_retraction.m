function [X_new, if_succeed] = conic_retraction(X, linear_sys, max_iter, pinf_eps, method)
    % retract X onto the spectrahedral cone {X psd : A X = 0}
    %   method = "alternating_projection" (default): Algorithm 3 of the paper
    %   method = "ssncg": orthogonal projection by semismooth Newton-CG
    if nargin < 5
        method = "alternating_projection";
    end
    switch method
        case "alternating_projection"
            [X_new, if_succeed] = conic_alternating_projection(X, linear_sys, max_iter, pinf_eps);
        case "ssncg"
            [X_new, if_succeed] = conic_ssncg_projection(X, linear_sys, max_iter, pinf_eps);
        otherwise
            error("Unknown retraction method: %s", method);
    end
end
