This folder contains the source code of the numerical tests described in the master thesis "Nonnegative matrix factorization and its application to Raman spectral data analysis" by Daniel Collin, Technische Universität Berlin. Find below a list of all files and instructions how to use the code. For further questions, contact the author: daniel.collin@web.de

----------------------------------------------------------

MATLAB files:
- ALS.m: function that computes an NMF with the ALS algorithm. Copyright: Lars Kai Hansen, IMM-DTU.
- beta_heuristic.m: function that computes a subset of frequencies based on the total activity over time.
- compute_K.m: function that computes the reaction coefficient matrix K.
- create_latex_plots.m: script that shows how to create some of the figures used in the thesis.
- eval_lorentz_sum.m: function that evaluates a sum of Lorentz functions.
- fit_lorentz.m: function that fits a sum of Lorentz functions to given data.
- gamma_heuristic.m: function that computes a subset of frequencies based on local maxima.
- generate_flow_network.m: function that generates a random flow network to be interpreted as a chemical reaction.
- generate_kinetics.m: function that generates a random reaction coefficient matrix K.
- generate_rand_params.m: function that generates random Lorentz parameters for component spectra W.
- generate_real_substance_data.m: function that generates a Raman spectral data set example, using real substances as species.
- get_H.m: function that generates the relative concentrations H based on K.
- greedy_lorentz.m: function that uses a greedy heuristic to approximate data by Lorentz functions.
- I.m: function that computes the interference of component spectra W.
- lorentz.m: function that returns a function handle to a Lorentz function.
- main.m: script that executes the numerical experiments described in the thesis.
- move_peaks.m: function that moves the base points of Lorentz functions closer together.
- NNLS.m: function that solves an NNLS problem with an exact block-coordinate descent scheme. Copyright: Nicolas Gillis, Université de Mons.
- NNLS_FPGM.m: function that solves an NNLS problem with a fast gradient method. Copyright: Nicolas Gillis, Université de Mons.
- NNLS_init.m: function that initializes the NNLS_FPGM solver. Copyright: Nicolas Gillis, Université de Mons.
- plot_mask.m: function that visualizes a selection of freqencies. 
- plot_results.m: function that plots selected results of the numerical experiment on the effects of increasing interference.
- raman_analysis.m: function that performs the analysis of a synthetic data set of time-resolved Raman spectral data.
- real_experiment.m: function that performs the analysis of real experimental data of time-resolved Raman spectral data.
- reorder_solution.m: function that permutes the NMF factors W and H such that their species are in the same order as the true factors.
- scale_factors.m: function that scales the NMF factors W and H such that they can be chemically interpreted.
- separability_experiment_1.m: function that executes the first numerical experiment on the effects of increasing interference.
- separability_experiment_2.m: function that executes the second numerical experiment on the effects of increasing interference.
- separability_experiment_3.m: function that executes the third numerical experiment on the effects of increasing interference
- SNPA.m: function that executes SNPA. Copyright: Nicolas Gillis, Université de Mons.
- SPA.m: function that executes SPA. Copyright: Nicolas Gillis, Université de Mons.
- XRAY.m: function that executes XRAY. Copyright: Nicolas Gillis, Université de Mons.


Other files:
- data: folder that contains data files of component spectra of real substances, the data set of real experimental data experimental_data.mat and the Lorentz parameters separable_p.mat used in the second experiment on the effect of increasing interferencs.
- README.txt: text file that contains a list of all files with short descriptions and instructions how to use the code.

----------------------------------------------------------

Instructions:

This code can be used to redo the numerical experiments described in the thesis and to reproduce the corresponding results, to generate the figures used in the thesis, and to do individually customized computations with the developed methodology.

To execute the numerical experiments, open the script main.m. Here, specify the desired settings for the experiments on the effect of increasing interference by using the options structure:

- options.n_exp: number of randomly generated component spectra (for experiment 1, default: 100) or reaction kinetics (for experiment 2, default: 100) or both (for experiment 3, default: 10).
- options.interference_level: levels of interference I(W) to be tested (default: I(W) = 0.0,0.1,...,1.0).
- options.noise_level: noise level delta that is used to disturb the data (default: delta = 0).
- options.algorithms: algorithms that are used to compute the NMF (1 = SPA, 2 = SNPA, 3 = XRAY, 4 = ALS, default: all for experiment 1 and 2, SNPA for experiment 3).
- options.n_species: numbers of species r for which the experiment is executed (only relevant for experiment 3, default: r = 2,3,4,5,6).

Then, by running the script once with the default settings for each experiment and once with the noise level delta = 0.01 for each experiment, the results shown in the thesis are reproduced. The results are saved in the folder results and then read and plotted by the function plot_results.

Additionally, one can do computations with one particular algorithm and data set, possessing a predefined interference and noise level. To do so, set for example options.n_exp = 1, options.interference_level = 0.25, options.noise_level = 0.005, options.algorithms = 3, options.n_species = 4 and execute separability_experiment_3. Then, one reaction of 4 species, consisting of random component spectra W with I(W) = 0.25 and reaction kinetics H, and disturbed by delta = 0.5% noise, is analyzed by XRAY. By opening the function raman_analysis.m and setting the variable PLOT = true, all steps of the analysis are visualized and stored in the folder figures.

The analysis of the real experimental data set can be done by running the second part of the script main.m. Here, specify a subinterval of frequencies if desired and set particular numbers of species r_range (default: r_range = 2,3,4,5,6) for which the analysis shall be performed. To visualize the results of the analysis and store the corresponding figures in the folder results, open the function real_experiment.m and set PLOT = true.

The script create_latex_plots.m shows how some of the figures used in the thesis can be generated. Execution leads to creation of the plots and its storage in the folder figures. This script can be modified as desired to create all sorts of plots.



