function mask = gamma_heuristic(A)
% This function selects a subset of frequencies by choosing only local
% maxima. It identifies almost linearly dependent rows close to intensity
% peaks.

[m,~] = size(A);

% Maximal value for each frequency.
M = max(A, [], 2);

% Identify local maxima in M.
ismax = zeros(m,1);
for k = 2:m-1
    if M(k-1) < M(k) && M(k+1) < M(k)
        ismax(k) = 1;
    end
end

mask = ismax;

end