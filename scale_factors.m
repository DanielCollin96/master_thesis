function [W, H] = scale_factors(W, H)
% This function scales the obtained NMF factors W and H such that the 
% product W*H remains unchanged, but now H is as close to a kinetics matrix
% (all row sums equal 1) as possible.

[r,n] = size(H);

% Set up linear program as decribed in chapter 4.1 to determine optimal 
% scaling factors d_i.
%   min vio_pos + vio_neg
%   s.t. H'*d + vio_pos - vio_neg = \ones
%        D*H <= 1 (entry wise)
%        ||W(:,j)||*d_i <= alpha*||W(:,i)||*d_j, for all i,j
%        d_i >= eps > 0
%        vio_pos, vio_neg >= 0

eps = 1e-9;
alpha = 5;

% Objective function.
f = [zeros(r,1); ones(n,1); ones(n,1)];

% Lower and upper bound for optimization variables.
lb = [ones(r,1) * eps; zeros(n,1); zeros(n,1)];
ub = inf * ones(r+2*n, 1);

% Equality conditions.
Aeq = [H', -speye(n), speye(n)];
beq = ones(n,1);

% Inequality conditions.
A1 = sparse((1:n*r)',kron(1:r,ones(1,n))',reshape(H',n*r,1),n*r,r+2*n);
A2 = spalloc(r^2-r,r+2*n,2*(r^2-r));
row_counter = 1;
for i = 1:(r-1)
    for j = (i+1):r
        factor1 = norm(W(:,j),'inf')/(norm(W(:,i),'inf')*alpha);
        factor2 = norm(W(:,i),'inf')/(norm(W(:,j),'inf')*alpha);
        A2(row_counter:row_counter+1,i) = [factor1; -1];
        A2(row_counter:row_counter+1,j) = [-1; factor2];
        row_counter = row_counter + 2;
    end
end
A = [A1; -A1; A2];
b = [ones(n*r, 1); zeros(n*r, 1); zeros(r^2-r,1)];

% Solve linear program.
options = optimoptions('linprog','ConstraintTolerance',1e-9,'Display','none');
[x,~,flag,~] = linprog(f,A,b,Aeq,beq,lb,ub,options);

% In case the conditions are too restrictive such that there is no feasible
% solution, mitigate conditions by increasing alpha.
while flag == -2
    alpha = alpha * 10;
    A2 = spalloc(r^2-r,r+2*n,2*(r^2-r));
    row_counter = 1;
    for i = 1:(r-1)
        for j = (i+1):r
            factor1 = norm(W(:,j),'inf')/(norm(W(:,i),'inf')*alpha);
            factor2 = norm(W(:,i),'inf')/(norm(W(:,j),'inf')*alpha);
            A2(row_counter:row_counter+1,i) = [factor1; -1];
            A2(row_counter:row_counter+1,j) = [-1; factor2];
            row_counter = row_counter + 2;
        end
    end
    A = [A1; -A1; A2];
    [x,~,flag,~] = linprog(f,A,b,Aeq,beq,lb,ub,options);
end
    
% Get scaling factors from solution.
d = x(1:r);
D = diag(d);
D_inv = diag(1.0./d);


% Perform scaling.
H = D*H;
W = W*D_inv;

end

