function interference = I(W)
% This function calculates the interference of component spectra W.

% Delete zero rows.
W(~any(W,2),:)= [];

[m,r] = size(W);

% Normalize rows.
W_norm = zeros(m,r);
for i = 1:m
    W_norm(i,:) = W(i,:)/norm(W(i,:),1);
end

% Search maximal entry of each column (each species).
val = zeros(r,1);
for i = 1:r
    val(i) = max(W_norm(:,i));
end

% Sum up deviations from perfect separability and scale to the interval
% [0,1].
interference = sum(ones(r,1)-val) / (r-1);

    
end