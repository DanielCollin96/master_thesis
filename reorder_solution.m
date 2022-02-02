function [W, H] = reorder_solution(W, H, W_true, H_true)
% This function reorders the rows and columns of the NMF factor H and W
% such that the ordering of species corresponds to the ordering of species
% of the true factor H_true and W_true.
% Note that this is a very weak heuristic and may be quite ineffective if 
% the given solution is not close to the true one.

% Ensure that the dimensions match and that all entries are numerical
% values.
assert(all(size(W) == size(W_true)));
assert(all(size(H) == size(H_true)));
assert(~ any(any(isnan(W))));
assert(~ any(any(isnan(H))));

r = size(W_true, 2);

% Compute similarity weights.
weights = W_true' * W;

% Initialize arrays for permutation.
perm = zeros(1,r);
zero_columns = all(weights == 0, 1);
non_used_idx = 1:r;

% Loop over all nonzero columns.
for s = 1:nnz(~zero_columns)
    
    % Find position of largest weight.
    [tmp, row_idx] = max(weights);
    [~, col_idx] = max(tmp);
    row_idx = row_idx(col_idx);
    
    % Store permutation.
    perm(col_idx) = row_idx;
    non_used_idx(row_idx) = 0;
    
    % Set used row and column to zero-
    weights(:, col_idx) = 0.0;
    weights(row_idx, :) = 0.0;
end

% Replace zero entries by non-used indices.
zero_entries = perm == 0;
non_used_idx = non_used_idx(non_used_idx ~= 0);
perm(zero_entries) = non_used_idx;

% Perform permutation of columns and rows of W and H such that they match
% to W_true and H_true.
W(:,perm) = W(:, 1:r);
H(perm, :) = H(1:r, :);

end