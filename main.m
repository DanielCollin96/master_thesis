% This script executes the numerical tests described in the master thesis
% "Nonnegative matrix factorization and its application to Raman spectral
% data analysis" by Daniel Collin, Technische Universität Berlin.
% Computed results and figures are stored in the folder "results".

if exist('results','dir') ~= 7
    mkdir('results');
end

%% Effect of increasing interference
rng(1)

% Initialize options for experiments.
options = struct;

% Number of randomly generated component spectra and/or reaction kinetics.
options.n_exp = 100;

% Noise level.
options.noise_level = 0;

% Levels of interference.
options.interference_level = 0:0.1:1; 

% Tested NMF algorithms. 1 = SPA, 2 = SNPA, 3 = XRAY, 4 = ALS.
options.algorithms = 1:4;

% Execute experiments and save results.

% Experiment 1 without noise.
result = separability_experiment_1(options);
save('results\separability_exp_1.mat','result')

% Experiment 1 with noise.
options.noise_level = 0.01;
result = separability_experiment_1(options);
save('results\separability_exp_1_noise.mat','result')

% Experiment 2 without noise.
options.noise_level = 0;
result = separability_experiment_2(options);
save('results\separability_exp_2.mat','result')

% Experiment 2 with noise.
options.noise_level = 0.01;
result = separability_experiment_2(options);
save('results\separability_exp_2_noise.mat','result')

% Experiment 3 without noise.
options.n_exp = 10;
options.noise_level = 0;
options.algorithms = 2;
result = separability_experiment_3(options);
save('results\separability_exp_3.mat','result')

% Experiment 3 with noise.
options.noise_level = 0.01;
result = separability_experiment_3(options);
save('results\separability_exp_3_noise.mat','result')


% Plot selected results.
plot_results(1,'separability_exp_1.mat')
plot_results(1,'separability_exp_1_noise.mat')

plot_results(2,'separability_exp_2.mat')
plot_results(2,'separability_exp_2_noise.mat')

plot_results(3,'separability_exp_3.mat')
plot_results(3,'separability_exp_3_noise.mat')



%% Real experiment

% Load experimental data.
data = load('data\experimental_data.mat');
A = data.A;
f = data.f;

% If desired, reduce to subinterval of frequencies.
% sub_ind = f >= 200 & f <= 750;
% f = f(sub_ind);
% A = A(sub_ind,:);

% Set interval of number of species that shall be computed.
r_range = 2:6;

% Execute data analysis.
[W_comp, H_comp, K_comp, error] = real_experiment(A,f,r_range)


