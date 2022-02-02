function error = raman_analysis(r, A, W_true, H_true, K_true, tspan, freq, noise_level, algorithms)
% This function performs the analysis of a data set of time-resolved Raman
% spectral data as described in chapter 4.1 and returns the relative 
% errors. In case there is specified a noise level, measurement noise is 
% added to the data set A and denoising is performed beforehand. To create 
% plots of the data and the results, set the variable PLOT = true. The 
% execution of the program is paused after every plot.

PLOT = false;

if exist('figures','dir') ~= 7
    mkdir('figures');
end

[m,n] = size(A);

% Add measurement noise and perform denoising if a positive noise level is
% specified.
if noise_level > 0
    
    % Add noise.
    A =  A + max(max(A)) * noise_level * abs(randn(size(A)));
    
    % Plot noisy data set.
    if PLOT
        surf(tspan,freq,A)
        axopts = {'FontSize', 18};
        fopts = {'Position', get(0, 'Screensize')};
        set(gca, axopts{:});
        set(gcf, fopts{:});
        ylabel('Wavenumber (cm^{-1})')%,'FontSize', 18);
        xlabel('Time (s)')%,'FontSize', 18);
        zlabel('Intensity')
        xlim([min(tspan) max(tspan)])
        ylim([min(freq) max(freq)])
        print(gcf,'figures\alg_example_a_noise.png','-dpng','-r300');
        pause
    end
    
    % Perform denoising by truncated SVD.
    [U,S,V] = svd(A);
    S_new = zeros(size(S));
    S_new(1:r,1:r) = S(1:r,1:r);
    A_trunc = U * S_new * V';
    A = max(A_trunc,0);
    
    % Plot denoised data set.
    if PLOT
        surf(tspan,freq,A)
        axopts = {'FontSize', 18};
        fopts = {'Position', get(0, 'Screensize')};
        set(gca, axopts{:});
        set(gcf, fopts{:});
        ylabel('Wavenumber (cm^{-1})')%,'FontSize', 18);
        xlabel('Time (s)')%,'FontSize', 18);
        zlabel('Intensity')
        xlim([min(tspan) max(tspan)])
        ylim([min(freq) max(freq)])
        print(gcf,'figures\alg_example_a_denoised.png','-dpng','-r300'); 
        pause
    end
else
    % If no noise is specified, plot original data set.
    if PLOT
        surf(tspan,freq,A)
        axopts = {'FontSize', 18};
        fopts = {'Position', get(0, 'Screensize')};
        set(gca, axopts{:});
        set(gcf, fopts{:});
        ylabel('Wavenumber (cm^{-1})')%,'FontSize', 18);
        xlabel('Time (s)')%,'FontSize', 18);
        zlabel('Intensity')
        xlim([min(tspan) max(tspan)])
        ylim([min(freq) max(freq)])
        print(gcf,'figures\alg_example_a.png','-dpng','-r300');
        pause
    end
end


% Selection of relevant frequencies.

% Beta heuristic, based on total activity over time.
beta = 2;
beta_mask = beta_heuristic(A, beta);
if PLOT
    plot_mask(A,beta_mask,freq,1)
    pause
end

% Gamma heuristic, based on local intensity maxima.
gamma_mask = gamma_heuristic(A);
if PLOT
    plot_mask(A,gamma_mask,freq,2)
    pause
end

% Combination of both masks.
mask = beta_mask & gamma_mask;
if PLOT
    plot_mask(A,mask,freq,3)
    pause
end

% If the heuristics eliminated too many frequencies, mitigate restrictions.
while nnz(mask) < r
    if nnz(gamma_mask) < r
        gamma_mask = true(m,1);
    else
        if beta > 1
            beta = beta * 0.5;
            beta_mask = beta_heuristic(A, beta);
        else
            beta_mask = true(m,1);
        end
    end
    mask = beta_mask & gamma_mask;
end

% Eliminate irrelevant frequencies.
selection = find(mask);
A_reduced = A(selection,:);

% Initialize structure to store the errors.
error = struct;

% Loop over all NMF algorithms. 
% 1 = SPA, 2 = SNPA, 3 = XRAY, 4 = ALS
for alg = algorithms
    
    % Compute near-separable NMF.
    if alg == 1
        subset_idx = SPA(A_reduced',r,1);
    elseif alg == 2
        options = struct;
        options.normalize = 1;
        options.display = 0;
        subset_idx = SNPA(A_reduced', r, options);
    elseif alg == 3
        subset_idx = XRAY(A_reduced', r);
    end
    
    % Compute H and W.
    if alg ~= 4
        
        % Transfer selected indices from reduced to complete data set.
        char_idx = selection(subset_idx);
        H = A(char_idx,:);
        
        % Solve NNLS problem to compute W.
        W = NNLS(A', H')';
        
        % Check if less than r columns have been extracted. If this is the
        % case, fill remaining rows/columns with zeros.
        r_comp = size(H,1);
        if r_comp < r
            H_comp = H;
            W_comp = W;
            H = zeros(r,n);
            W = zeros(m,r);
            H(1:r_comp,:) = H_comp;
            W(:,1:r_comp) = W_comp;
        end
        
        % Define the first computed factors as initial matrices for ALS.
        if exist('W_init','var') == 0
            W_init = W;
            H_init = H;
        end
    elseif alg == 4
        
        % If no other factors have been computed previously, initialize ALS
        % with random matrices.
        if exist('W_init','var') == 0
            W_init = rand(m,r);
            H_init = rand(r,n);
        end
        
        % Compute NMF with ALS.
        [W,H] = ALS(A, r, 500, 0, W_init, H_init);
        
        % Keep solution only if ALS improves the approximation.
        error_pre = norm(A - W_init*H_init,'fro')/norm(A,'fro');
        error_post = norm(A - W*H,'fro')/norm(A,'fro');
        if error_pre < error_post
            W = W_init;
            H = H_init;
        end
    end

    % Get zero rows and columns of the factors to exclude them from the
    % following scaling.
    zero_rows = all(H == 0, 2);
    zero_columns = all(W == 0, 1)';
    zero = zero_rows & zero_columns;
    r_comp = nnz(~zero);
    
    % Scale NMF factors such that H is close a kinetics matrix.
    [W(:,~zero),H(~zero,:)] = scale_factors(W(:,~zero), H(~zero,:));
    
    % Reorder rows and columns of H and W such that they have the same
    % order as the rows and columns of the true matrices. This is necessary
    % to compute the relative errors.
    [W,H] = reorder_solution(W, H, W_true, H_true);
    assert( ~any(isnan(W(:))) && ~any(isinf(W(:))) );
     
    if PLOT
        % Plot comparison of computed and true component spectra as well as
        % the characteristic frequencies.
        for i = 1:r
            plot(freq,W(:,i),'LineWidth',0.6)
            hold on
            plot(freq,W_true(:,i),'k--','LineWidth',0.6)
            yl = ylim;
            for j = 1:r
                line([freq(char_idx(j)) freq(char_idx(j))], [yl(1),yl(2)],'LineWidth',0.6, 'LineStyle', '-', 'Color', [0, 0, 0, 0.4]);
            end
            hold off
            axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
            fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
            set(gca, axopts{:});
            set(gcf, fopts{:});
            xlabel('Wavenumber (cm^{-1})')%,'FontSize', 18);
            ylabel('Intensity')%,'FontSize', 18);
            xlim([min(freq) max(freq)])
            legend({'Computed','Original'},'Location','northwest')
            print(gcf,strcat('figures\alg_example_w',int2str(i),'.png'),'-dpng','-r300'); 
            pause
        end
    
        % Plot comparison of computed and true reaction kinetics.
        plot(tspan,H','LineWidth',0.6)
        hold on
        plot(tspan,H_true','k--','LineWidth',0.6)
        hold off
        axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
        fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
        set(gca, axopts{:});
        set(gcf, fopts{:});
        xlabel('Time (s)')%,'FontSize', 18);
        ylabel('Relative concentration')%,'FontSize', 18);
        xlim([min(tspan) max(tspan)])
        legend_names = {'Species A','Species B','Species C','Species D','Species E','Species F','Species G','Species H'};
        legend_names = {legend_names{1:r},'Original'};
        legend(legend_names,'Location','northwest')
        print(gcf,'figures\alg_example_h.png','-dpng','-r300'); 
        pause
    end
    
    % Compute reaction coefficients K.
    K = compute_K(H,tspan);

    % Plot comparison of the kinetics based on the computed reaction
    % coefficients and the true kinetics.
    if PLOT
        H_K = get_H(K,tspan,H(:,1));
        
        plot(tspan,H_K','LineWidth',0.6)
        hold on
        plot(tspan,H_true','k--','LineWidth',0.6)
        hold off
        axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
        fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
        set(gca, axopts{:});
        set(gcf, fopts{:});
        xlabel('Time (s)')%,'FontSize', 18);
        ylabel('Relative concentration')%,'FontSize', 18);
        xlim([min(tspan) max(tspan)])
        legend_names = {'Species A','Species B','Species C','Species D','Species E','Species F','Species G','Species H'};
        legend_names = {legend_names{1:r},'Original'};
        legend(legend_names,'Location','northwest')
        print(gcf,'figures\alg_example_k.png','-dpng','-r300'); 
        pause
    end
    
    % Compute relative errors and store them in error structure.
    result = zeros(5,1);
    result(1) = norm(A - W*H,'fro')/norm(A,'fro');
    result(2) = norm(W_true - W,'fro')/norm(W_true,'fro');
    result(3) = norm(H_true - H,'fro')/norm(H_true,'fro');
    result(4) = norm(K_true - K,'fro')/norm(K_true,'fro');
    result(5) = r_comp;
    
    if alg == 1
        error.SPA = result;
    elseif alg == 2
        error.SNPA = result;       
    elseif alg == 3
        error.XRAY = result;
    elseif alg == 4
        error.ALS = result;
    end
    
end
