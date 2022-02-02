function mask = beta_heuristic(A, beta)
% This function selects a subset of frequencies based on the total activity
% over time. Only those frequencies are selected, whose activity is at 
% least beta * (minimum activity).

% Ensure beta is at least 1.
if beta < 1
    beta = 1;
end

m = size(A,1);

% Determine minimum activity over all frequencies.
rowsum = sum(abs(A), 2);
min_activity = min(rowsum);

% Only select frequencies with sufficient activity.
mask = false(m,1);
mask(rowsum > beta * min_activity) = true;


end