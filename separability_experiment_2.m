function result = separability_experiment_2(options)
% This function executes the second numerical experiment with randomly
% generated data. It computes the relative approximation errors of A, W, H
% and K. Here, the same initial compononent spectra are used for every data
% set and n_exp reaction kinetics are generated.

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

% Number of species.
n_species = 5;

% Range of discrete frequencies.
f = linspace(500,3000,2501);

% Range of discrete time points.
tspan = linspace(0,50,151);

% Load Lorentz parameters for almost separable component spectra.
data = load('data\separable_p.mat');
p_orig = data.p_orig;

% Create corresponding W component spectra.
W_orig = zeros(length(f),n_species);
for i = 1:n_species
    W_orig(:,i) = eval_lorentz_sum(p_orig{i},f);
end

% Compute interference of original component spectra.
interference_orig = I(W_orig);

% Randomly generate reaction kinetics.
K_storage = cell(n_exp,1);
H_storage = cell(n_exp,1);
for i = 1:n_exp
    
    % Generate random K and h0.
    [K,h0] = generate_kinetics(n_species,0.2);
    
    % Compute corresponding H.
    H = get_H(K,tspan,h0);
    
    % Compute singular values and store matrices only if the smallest
    % singular value is greater than 0.2 to avoid bad conditioned matrices.
    % In case of bad condition, compute new kinetics.
    S = svd(H);
    while min(S) < 0.2
        [K,h0] = generate_kinetics(n_species,0.2);
        H = get_H(K,tspan,h0);
        S = svd(H);
    end
    K_storage{i} = K;
    H_storage{i} = H;
end


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

% Accuracy treshold for bisection method.
eps = 1e-6;


% Loop over all interference levels.
for k = 1:length(interference_level)
    
    fprintf('Computing NMF for separability = %.2f.\n',interference_level(k))
   
    % Initial interference and component spectra.
    current_interference = interference_orig;
    W_true = W_orig;

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

    % Store actual value of I(W) and lambda. Note that these values are
    % the same for each data set of one interference level because the
    % same component spectra are used.
    average_interference(k) = current_interference;
    average_lambda(k) = lam;
        
    % Loop over all data sets for each level of interference.
    for i = 1:n_exp
        interference(counter) = current_interference;
        
        % Compute measurement matrix.
        H_true = H_storage{i};
        K_true = K_storage{i};
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
    
    fprintf('Lambda: %.4f.\n',lambda(k));
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