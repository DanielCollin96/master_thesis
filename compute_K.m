function K = compute_K(H,tspan)
% This function computes the reaction coefficient matrix K, given the
% computed relative concentrations H, by solving the linear least squares
% problem 
%       min ||\mathcal{D}-\tilde{\mathcal{H}}\tilde{\kappa}||_F^2,
% as described in chapter 4.1.

[r,n] = size(H);

% Set up matrix \tilde{\mathcal{H}}. Exclude the variables k_{ss} and the
% columns belonging to it. Insert the vectors −H(i,:)'.
C = spalloc(r*(n-1),r*(r-1),2*r*(r-1)*(n-1));
for i = 1:r
    row_range = (i-1)*(n-1)+1:i*(n-1);
    col_range = (i-1)*(r-1)+1:i*(r-1);
    
    % Exclude i-th column and insert all other columns at their
    % corresponding place.
    H_range = [1:i-1, i+1:r];
    C(row_range,col_range) = H(H_range,1:n-1)';
    for j = 1:r
        if i ~= j
            col = (j-1)*(r-1);
            if i < j
                col = col + i;
            else % i > j
                col = col + i-1;
            end
            C(row_range,col) = - H(i,1:n-1)';
        end
    end
end

% Set bounds for reaction coefficients.
lb = zeros(r*(r-1),1);
ub = ones(r*(r-1),1);

% Compute \mathcal{D} as difference quotient with interpolated values.
h = 1e-8;
H_interp = zeros(r,n-1);
for i = 1:r
    H_interp(i,:) = spline(tspan,H(i,:),tspan(1:n-1)+h);
end
H_dif = H_interp - H(:,1:n-1);

d = reshape(H_dif',r*(n-1),1);
d = d/h;

% Solve linear least squares problem.
options = optimoptions('lsqlin','Display','off');
[K,~,~,~,~,~] = lsqlin(C,d,[],[],[],[],lb,ub,[],options);

% Bring solution (currently in vector form) into matrix form and compute
% diagonal entries.
K_temp = reshape(K,r-1,r)';
L = [tril(K_temp,-1), zeros(r,1)];
U = [zeros(r,1), triu(K_temp)];
K = L + U;
for i = 1:r
    K(i,i) = -sum(K(:,i));
end     

end