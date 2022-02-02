function p = move_peaks(p,freq,n_meet,lambda)
% This function moves the base points of Lorentz functions closer together
% by convex combinations. A number of n_meet focal points is distributed
% equidistantly over the interval of frequencies and each base point is
% moved towards the closest focal point.

% Determine distance between the focal points.
freq_min = min(freq);
freq_max = max(freq);
interval_length = (freq_max - freq_min)/(n_meet+1);

% Calculate focal points.
meet_points = zeros(n_meet,1);
meet_points(1) = freq_min + interval_length;
for i = 2:n_meet
    meet_points(i) = meet_points(i-1) + interval_length;
end

% Loop over all species.
n_species = size(p,1);
for i = 1:n_species
    
    % Loop over all Lorentzians.
    n_lor = size(p{i},1);
    for j = 1:n_lor
        
        % Get current base point.
        x0 = p{i}(j,1);
        
        % Determine closest focal point.
        mp = round((x0 - freq_min) / interval_length);
        if mp == 0
            mp = 1;
        elseif mp == n_meet + 1
            mp = n_meet;
        end
        closest_point = meet_points(mp);
        
        % Compute convex combination.
        x0 = (1-lambda) * x0 + lambda * closest_point;
        p{i}(j,1) = x0;
    end
end

end