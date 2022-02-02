function y = eval_lorentz_sum(params, x)
% This function evaluates a sum of Lorentzians defined through the 
% parameter array params at x.
% Each row of 'params' has length 3 with the following meaning:
% - params(:,1): x0, the base point
% - params(:,2): gamma, the half-width at half height
% - params(:,3): I, the intensity


[n_lor, n] = size(params);
assert(n==3);

y = zeros(size(x));

% Loop over all Lorentz functions and add values at x.
for k=1:n_lor
    x0 = params(k,1);
    gamma = params(k,2);
    I = params(k,3);

    f = lorentz(x0,gamma,I);
    y = y + f(x);
end

end