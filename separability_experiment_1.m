function result = separability_experiment_1(options)
% This function executes the first numerical experiment with randomly
% generated data. It computes the relative approximation errors of A, W, H
% and K. Here, the same reaction kinetics is used for every data set and 
% n_exp component spectra are generated.

% Set options.
if isfield(options,'n_exp')
    n_exp = options.n_exp;
else
    n_exp = 100;
end

if isfield(options,'interference_level')
    interference_level = options.interference_level;
else
    interference_level = 0:0.1:1;
end

if isfield(options,'noise_level')
    noise_level = options.noise_level;
else
    noise_level = 0;
end

if isfield(options,'algorithms')
    algorithms = options.algorithms;
else
    % 1 = SPA, 2 = SNPA, 3 = XRAY, 4 = ALS.
    algorithms = 1:4;
end

% Number of focal points used for moving intensity peaks of spectra closer 
% together.
n_meet = 1;

% Reaction kinetics matrix K.
K_true = [...
    -0.53,  0.02,  0.0,   0,  0.0;
     0.53, -0.66,  0.25,  0,  0.0;
     0.0,   0.43, -0.36,  0,  0.1;
     0.0,   0.21,  0.0,   0,  0.0;
     0.0,   0.0,   0.11,  0, -0.1];

% Number of species.
n_species = size(K_true,1);

% Initial relative concentrations.
h0 = zeros(n_species,1);
h0(1) = 1;

% Range of discrete frequencies.
f = linspace(500,3000,2501);

% Range of discrete time points.
tspan = linspace(0,50,151);

% Compute relative concentrations H based on K.
H_true = get_H(K_true,tspan,h0);

% Initialize arrays to store results.

% Value of interference I(W) for component spectra.
interference = zeros(1,length(interference_level)*n_exp);

% Value of interference I(W) for component spectra, averaged over all W for
% each level of interference.
average_interference = zeros(1,length(interference_level));

% Value of lambda of the convex combination used for moving the intensity
% peaks of a spectrum closer together.
lambda = zeros(1,length(interference_level)*n_exp);

% Value of lambda of the convex combination used for moving the intensity
% peaks of a spectrum closer together, averaged over all W for each
% level of interference.
average_lambda = zeros(1,length(interference_level));

% Various results for each algorithm separately. Meaning of the rows:
% 1: relative approximation error of A
% 2: relative approximation error of W
% 3: relative approximation error of H
% 4: relative approximation error of K
% 5: number of computed species
single = struct;
single.SPA = zeros(5,length(interference_level)*n_exp);
single.SNPA = zeros(5,length(interference_level)*n_exp);
single.XRAY = zeros(5,length(interference_level)*n_exp);
single.ALS = zeros(5,length(interference_level)*n_exp);

% Various results for each algorithm separately, averaged for each level of
% interference.
average = struct;
average.SPA = zeros(5,length(interference_level));
average.SNPA = zeros(5,length(interference_level));
average.XRAY = zeros(5,length(interference_level));
average.ALS = zeros(5,length(interference_level));

% Counter pointing to entry where the next result is stored.
counter = 1;

% Array to store randomly generated Lorentz parameters.
params = cell(n_exp,1);

% Generate random Lorentz parameters for component spectra.
for i = 1:n_exp
    params{i} = generate_rand_params(f,n_species);
end

% Accuracy treshold for bisection method.
eps = 1e-6;

% Loop over all interference levels.
for k = 1:length(interference_level)
    
    fprintf('Computing NMF for separability = %.2f.\n',interference_level(k))
    
    % Loop over all data sets for each level of interference.
    for i = 1:n_exp
        
        p_orig = params{i};
        
        % Create W for initial Lorentz parameters.
        W_true = zeros(length(f),n_species);
        for j = 1:n_species
            W_true(:,j) = eval_lorentz_sum(p_orig{j},f);
        end
        
        % Measure interference.
        current_interference = I(W_true);
        
        % Set initial upper bound, lower bound and lambda for bisection 
        % method. 
        ub = 1;
        lb = 0;
        lam = 0;
        
        % Count iterations.
        iteration = 0;
        
        % Bisection method to generate W with specific I(W). Limit to 100
        % iterations.
        while abs(current_interference - interference_level(k)) > eps && iteration <= 100
            
            % Bisect interval depending on I(W) of the current parameters.
            if current_interference > interference_level(k)
                ub = lam;
            else
                lb = lam;
            end
            
            % Set lambda to the middle of the interval.
            lam = (ub + lb)/2;
            
            % Move intensity peaks by convex combination with new lambda.
            p = move_peaks(p_orig,f,n_meet,lam);
            
            % Create new W for modified Lorentz parameters.
            for j = 1:n_species
                W_true(:,j) = eval_lorentz_sum(p{j},f);
            end
            
            % I(W) of modified component spectra.
            current_interference = I(W_true);
            
            iteration = iteration + 1;
        end
        
        % Store actual value of I(W) and lambda.
        interference(counter) = current_interference;
        lambda(counter) = lam;
        
        % Compute measurement matrix.
        A = W_true * H_true;
        
        % Execute the computational analysis of the data set and store the
        % errors.
        error = raman_analysis(n_species, A, W_true, H_true, K_true, tspan, f, noise_level, algorithms);

        if isfield(error,'SPA')
            single.SPA(:,counter) = error.SPA;
        end
        if isfield(error,'SNPA')
            single.SNPA(:,counter) = error.SNPA;
        end
        if isfield(error,'XRAY')
            single.XRAY(:,counter) = error.XRAY;
        end
        if isfield(error,'ALS')
            single.ALS(:,counter) = error.ALS;
        end
        
        counter = counter + 1;
    end
    
    % Compute averages of results for each level of interference.
    average_interference(k) = sum(interference((k-1)*n_exp+1:n_exp*k)) / n_exp;
    average_lambda(k) = sum(lambda((k-1)*n_exp+1:n_exp*k)) / n_exp;
    
    if isfield(error,'SPA')
        average.SPA(:,k) = sum(single.SPA(:,(k-1)*n_exp+1:n_exp*k),2) / n_exp;
    end
    if isfield(error,'SNPA')
        average.SNPA(:,k) = sum(single.SNPA(:,(k-1)*n_exp+1:n_exp*k),2) / n_exp;
    end
    if isfield(error,'XRAY')
        average.XRAY(:,k) = sum(single.XRAY(:,(k-1)*n_exp+1:n_exp*k),2) / n_exp;
    end
    if isfield(error,'ALS')
        average.ALS(:,k) = sum(single.ALS(:,(k-1)*n_exp+1:n_exp*k),2) / n_exp;
    end
    
    fprintf('Average lambda: %.4f.\n',average_lambda(k));
end

% Store results in structure.
result = struct;
result.single = single;
result.average = average;
result.interference = interference;
result.average_interference = average_interference;
result.lambda = lambda;
result.average_lam = average_lambda;
result.noise = noise_level;

end