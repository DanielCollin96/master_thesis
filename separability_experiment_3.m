function result_array = separability_experiment_3(options)
% This function executes the second numerical experiment with randomly
% generated data. It computes the relative approximation errors of A, W, H
% and K. Here, the same initial compononent spectra are used for every data
% set and n_exp reaction kinetics are generated.

% Set options.
if isfield(options,'n_exp')
    n_exp = options.n_exp;
else
    n_exp = 10;
end

if isfield(options,'n_species')
    n_species = options.n_species;
else
    n_species = 2:6;
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
    algorithms = 2;
end

% Number of focal points used for moving intensity peaks of spectra closer 
% together.
n_meet = 1;

% Range of discrete frequencies.
f = linspace(500,3000,2501);

% Range of discrete time points.
tspan = linspace(0,50,151);

% Generate random Lorentz parameters for each number of species.
params = cell(length(n_species),n_exp);
for i = 1:length(n_species)
    for j = 1:n_exp
        params{i,j} = generate_rand_params(f,n_species(i));
    end
end

% Generate random reaction kinetics for each number of species.
K_storage = cell(length(n_species),n_exp);
H_storage = cell(length(n_species),n_exp);
for i = 1:length(n_species)
    for j = 1:n_exp
        % Generate random K and h0.
        [K,h0] = generate_kinetics(n_species(i),0.2);

        % Compute corresponding H.
        H = get_H(K,tspan,h0);

        % Compute singular values and store matrices only if the smallest
        % singular value is greater than 0.2 to avoid bad conditioned matrices.
        % In case of bad condition, compute new kinetics.
        S = svd(H);
        while min(S) < 0.2
            [K,h0] = generate_kinetics(n_species(i),0.2);
            H = get_H(K,tspan,h0);
            S = svd(H);
        end
        
        K_storage{i,j} = K;
        H_storage{i,j} = H;
    end
end

% Accuracy treshold for bisection method.
eps = 1e-6;

% Initialize array to store results.
result_array = cell(length(n_species),1);

% Loop over all numbers of species.
for i = 1:length(n_species)
    
    % Set current number of species.
    r = n_species(i);
    
    fprintf('Start computations for r = %i.\n',r);
    
    % Initialize arrays to store results for each number of species.

    % Value of interference I(W) for component spectra.
    interference = zeros(1,length(interference_level)*n_exp^2);

    % Value of interference I(W) for component spectra, averaged over all W for
    % each level of interference.
    average_interference = zeros(1,length(interference_level));

    % Value of lambda of the convex combination used for moving the intensity
    % peaks of a spectrum closer together.
    lambda = zeros(1,length(interference_level)*n_exp^2);

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
    single.SPA = zeros(5,length(interference_level)*n_exp^2);
    single.SNPA = zeros(5,length(interference_level)*n_exp^2);
    single.XRAY = zeros(5,length(interference_level)*n_exp^2);
    single.ALS = zeros(5,length(interference_level)*n_exp^2);

    % Various results for each algorithm separately, averaged for each level of
    % interference.
    average = struct;
    average.SPA = zeros(5,length(interference_level));
    average.SNPA = zeros(5,length(interference_level));
    average.XRAY = zeros(5,length(interference_level));
    average.ALS = zeros(5,length(interference_level));
    
    % Counter pointing to entry where the next result is stored.
    counter = 1;
    
    % Loop over all levels of interference.
    for j = 1:length(interference_level)
        
        % Loop over all spectra.
        for spec = 1:n_exp
            
            % Create W for initial Lorentz parameters.
            p_orig = params{i,spec};
            W_true = zeros(length(f),r);
            for col = 1:r
                W_true(:,col) = eval_lorentz_sum(p_orig{col},f);
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
            while abs(current_interference - interference_level(j)) > eps && iteration <= 100

                % Bisect interval depending on I(W) of the current parameters.
                if current_interference > interference_level(j)
                    ub = lam;
                else
                    lb = lam;
                end

                % Set lambda to the middle of the interval.
                lam = (ub + lb)/2;

                % Move intensity peaks by convex combination with new lambda.
                p = move_peaks(p_orig,f,n_meet,lam);
                
                % Create new W for modified Lorentz parameters.
                for col = 1:r
                    W_true(:,col) = eval_lorentz_sum(p{col},f);
                end
                
                % I(W) of modified component spectra.
                current_interference = I(W_true);
                
                iteration = iteration + 1;
            end
            
            % Loop over all kinetics.
            for kin = 1:n_exp
                % Store actual value of I(W) and lambda.
                interference(counter) = current_interference;
                lambda(counter) = lam;
                
                % Compute measurement matrix.
                H_true = H_storage{i,kin};
                K_true = K_storage{i,kin};
                A = W_true * H_true;
                
                % Execute the computational analysis of the data set and store the
                % errors.
                error = raman_analysis(r, A, W_true, H_true, K_true, tspan, f, noise_level, algorithms);
        
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
        end
        
        % Compute averages of results for each level of interference.
        average_interference(j) = sum(interference((j-1)*n_exp^2+1:n_exp^2*j)) / n_exp^2;
        average_lambda(j) = sum(lambda((j-1)*n_exp^2+1:n_exp^2*j)) / n_exp^2;
        
        if isfield(error,'SPA')
            average.SPA(:,j) = sum(single.SPA(:,(j-1)*n_exp^2+1:n_exp^2*j),2) / n_exp^2;
        end
        if isfield(error,'SNPA')
            average.SNPA(:,j) = sum(single.SNPA(:,(j-1)*n_exp^2+1:n_exp^2*j),2) / n_exp^2;
        end
        if isfield(error,'XRAY')
            average.XRAY(:,j) = sum(single.XRAY(:,(j-1)*n_exp^2+1:n_exp^2*j),2) / n_exp^2;
        end
        if isfield(error,'ALS')
            average.ALS(:,j) = sum(single.ALS(:,(j-1)*n_exp^2+1:n_exp^2*j),2) / n_exp^2;
        end
    end
    
    result = struct;
    result.single = single;
    result.average = average;
    result.interference = interference;
    result.average_interference = average_interference;
    result.lambda = lambda;
    result.average_lambda = average_lambda;
    result.noise = noise_level;
    result.r = r;
    result_array{i} = result;
end

end