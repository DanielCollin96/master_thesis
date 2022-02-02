function p = generate_rand_params(freq,n_species)
% This function generates random Lorentz parameters for component spectra.

p = cell(n_species,1);

freq_min = min(freq);
freq_max = max(freq);

for i = 1:n_species
    
    % Determine number of Lorentzians the spectrum consists of.
    n_lor = randi([5,15]);
    p{i} = zeros(n_lor,3);
    
    % Determine parameters for each function.
    for j = 1:n_lor
        x0 = randi([freq_min,freq_max]);
        sigma = 1 + abs(randn) * 5;     
        I = 100 + abs(randn) * 10000;
        p{i}(j,:) = [x0,sigma,I];
    end
end

end