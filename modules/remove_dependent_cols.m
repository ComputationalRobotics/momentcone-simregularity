function [Aind, ids] = remove_dependent_cols(A, tau)
    if nargin < 2
        tau = max(sqrt(eps), size(A,2)*eps);   % default tolerance
    end
    if size(A, 2) < 2
        Aind = A;
        ids = 1;
        return;
    end
    % QR with column pivoting
    [~, R, p] = qr(A, 0);                      % economy size
    % find effective rank
    diagR = abs(diag(R));
    r = find(diagR > tau*diagR(1), 1, 'last');
    if isempty(r), r = 0; end
    % extract independent columns
    ids = p(1:r);
    Aind = A(:, ids);
end
