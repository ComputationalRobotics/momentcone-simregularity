clear; close all;
rng(0);

%% add paths
pathinfo = dictionary();

pathinfo("mosek") = "~/ksc/matlab-install/mosek/10.1/toolbox/r2017a";
pathinfo("msspoly") = "~/ksc/matlab-install/spotless";
pathinfo("sdpt3") = "~/ksc/matlab-install/SDPT3-4.0";

pathinfo("sparsesdprelax") = "./sos-sdp-conversion";
pathinfo("modules") = "./modules";

keys = pathinfo.keys;
for i = 1: length(keys)
    key = keys(i);
    addpath(genpath(pathinfo(key)));
end

%% hyper-parameters
n = 9;
kappa = 2;
s = 34; % z_i's number
p = 101; % prime number for sample 

%% special case: kappa = 2
if kappa == 2
    tmp1 = nchoosek(n+2, 2);
    tmp2 = nchoosek(n+4, 4);
    A_poly_row_num_theory = 0.5 * tmp1 * (tmp1 + 1) - tmp2;
    a = A_poly_row_num_theory;
    s_max = floor(0.5 * (1 + sqrt(1 + 8*a)));
end

%% set Mosek parameters
param = struct(); 
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-9;  % objective gap
param.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-15;  % primal feasibility
param.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-13;  % dual   feasibility
param.MSK_DPAR_INTPNT_CO_TOL_INFEAS  = 1e-15; % infeasibility test

%% generate moment constraints
mat_size = nchoosek(n + kappa, kappa);
[At_sdpt3, others] = generate_moment_cone(n, kappa, false);
At_sedumi = others.At_sedumi;
m = size(At_sedumi, 2);
% force moment constraints to be all integers
At_mom = round(2 * At_sedumi);
% for i = 1: size(At_mom, 2)
%     Atmp = At_mom(:, i);
%     Atmp = reshape(Atmp, mat_size, mat_size);
% end

%% construnct Z and V     
% {z_1, ..., z_s} is of size n by 1
% entries iid uniform in {0,...,p-1}
Z = randi([0, p-1], n, s);   
x = msspoly('x', n);
v_sym = mymonomial(x, kappa);
V = zeros(size(v_sym, 1), s);
for i = 1: s
    V(:, i) = double(subs(v_sym, x, Z(:, i)));
end
% check the full column rankness over GF_p
% r = gfrank(mod(V, p), p);
r = rank_modp(V, p);
isFullColRank = (r == size(V,2));
if ~isFullColRank 
    fprintf("V itself is not full column rank over GF(%d)! \n", p);
    return;
end

%% generate polyhedral constraints A_poly
A_poly_column_num = round(0.5 * s * (s-1));
A_poly_row_num = size(At_mom, 2);
A_poly = zeros(A_poly_row_num, A_poly_column_num);
for i = 1: size(At_mom, 2)
    Atmp = At_mom(:, i); 
    Atmp = reshape(Atmp, mat_size, mat_size);
    tmp = Atmp * V; tmp = V' * tmp;
    tmp = tmp(triu(true(size(tmp)), 1)).';
    A_poly(i, :) = tmp;
end
% check column rank of A_poly 
r = rank_modp(A_poly, p);
isFullColRank = (r == size(A_poly,2));
fprintf("n = %d, kappa = %d, s = %d \n", n, kappa, s);
if isFullColRank
    fprintf("column rank over GF(%d) = column number = %d \n", p, r);
else
    fprintf("A_poly is not full column rank over GF(%d)! \n", p);
    fprintf("column rank over GF(%d) is %d, but column number is %d ... \n", p, r, size(A_poly,2));
    return;
end


%% helper function
function rk = rank_modp(A, p)
    % Exact rank over GF(p) for prime p
    A = mod(double(A), p);
    [m,n] = size(A);
    rk = 0; r = 1;
    
    for c = 1:n
        % find pivot at/below row r
        pivot = find(A(r:m,c)~=0, 1, 'first');
        if isempty(pivot), continue; end
        pivot = pivot + r - 1;
    
        % swap pivot row up
        if pivot ~= r
            A([r pivot],:) = A([pivot r],:);
        end
    
        % normalize pivot row
        invPivot = powermod(A(r,c), p-2, p);   % Fermat inverse since p is prime
        A(r,:) = mod(A(r,:) * invPivot, p);
    
        % eliminate column c in all other rows
        for i = 1:m
            if i ~= r && A(i,c) ~= 0
                A(i,:) = mod(A(i,:) - A(i,c)*A(r,:), p);
            end
        end
    
        rk = rk + 1;
        r = r + 1;
        if r > m, break; end
        if rk == n, break; end
    end
end

function y = powermod(a,e,p)
    y = 1; a = mod(a,p);
    while e > 0
        if bitand(e,1), y = mod(y*a,p); end
        a = mod(a*a,p);
        e = bitshift(e,-1);
    end
end



