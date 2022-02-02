function [W_comp, H_comp, K_comp, error] = real_experiment(A, f, r_range)
% This function performs the analysis of real experimental data of 
% time-resolved Raman spectral data as described in chapter 4.3 for a range
% of different numbers of species and returns the computed factors and the
% relative approximation error of the data. To create plots of the data and
% the results, set the variable PLOT = true.

PLOT = true;

if exist('results\real_experiment','dir') ~= 7
    mkdir('results\real_experiment');
end


[m,n] = size(A);
tspan = linspace(0,42*(n-1),n);

% Plot original data set.
if PLOT
    surf(tspan,f,A)
    axopts = {'FontSize', 18};
    fopts = {'Position', get(0, 'Screensize')};
    set(gca, axopts{:});
    set(gcf, fopts{:});
    ylabel('Wavenumber (cm^{-1})')%,'FontSize', 18);
    xlabel('Time (s)')%,'FontSize', 18);
    zlabel('Intensity')
    xlim([min(tspan) max(tspan)])
    ylim([min(f) max(f)])
    print(gcf,'results\real_experiment\a_orig.png','-dpng','-r300'); 
    pause
end

% Compute buffer spectrum.
buffer = min(A, [], 2);

% Plot buffer spectrum.
if PLOT
    plot(f,buffer,'LineWidth',0.6)
    axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
    fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
    set(gca, axopts{:});
    set(gcf, fopts{:});
    xlabel('Wavenumber (cm^{-1})')%,'FontSize', 18);
    ylabel('Intensity')%,'FontSize', 18);
    xlim([min(f) max(f)])
    %ylim([0 max(max(W))*1.01])
    print(gcf,'results\real_experiment\buffer.png','-dpng','-r300'); 
    pause
end

% Subtract buffer spectrum from data.
offset = repmat(buffer, 1, n);
A = A - offset;
A_without_buffer = A;


% Initialize structures for results.
num_r = length(r_range);
error = zeros(num_r,1);
W_comp = cell(num_r,1);
H_comp = cell(num_r,1);
K_comp = cell(num_r,1);

% Analyze data set for each number of species.
counter = 1;
for r = r_range
    
    % Perform denoising by truncated SVD.
    [U,S,V] = svd(A_without_buffer);
    S_new = zeros(size(S));
    S_new(1:r,1:r) = S(1:r,1:r);
    A_trunc = U * S_new * V';
    A = max(A_trunc,0);
    
    if PLOT
        surf(tspan,f,A)
        axopts = {'FontSize', 18};
        fopts = {'Position', get(0, 'Screensize')};
        set(gca, axopts{:});
        set(gcf, fopts{:});
        ylabel('Wavenumber (cm^{-1})')%,'FontSize', 18);
        xlabel('Time (s)')%,'FontSize', 18);
        zlabel('Intensity')
        xlim([min(tspan) max(tspan)])
        ylim([min(f) max(f)])
        print(gcf,'results\real_experiment\a_preprocessed.png','-dpng','-r300'); 
        pause
    end
    
    % Selection of relevant frequencies.
    beta_mask = beta_heuristic(A, 5);
    gamma_mask = gamma_heuristic(A);
    mask = beta_mask & gamma_mask;
    fprintf('Selected mask covers %d frequencies.\n', nnz(mask));
   
    % Eliminate irrelevant frequencies.
    selection = find(mask);
    A_reduced = A(selection, :);
    
    % Compute near-separable NMF.
    options = struct;
    options.normalize = 1;
    options.display = 0;
    subset_idx = SNPA(A_reduced', r, options);

    % Transfer selected indices from reduced to complete data set.
    char_idx = selection(subset_idx);
    H = A(char_idx,:);
    
    % Solve NNLS problem to compute W.
    W = NNLS(A', H')';

    % Compute NMF with ALS.
    [W_new,H_new] = ALS(A, 5, 500, 0, W, H);
    
    % Keep solution only if ALS improves the approximation.
    error_pre = norm(A - W*H,'fro') / norm(A,'fro');
    error_post = norm(A - W_new*H_new,'fro') / norm(A,'fro');
    if error_post < error_pre
        W = W_new;
        H = H_new;
        disp('Local improvement, using ALS solution.');
    else
        disp('No local improvement, using SNPA solution.');
    end
    
    % Scale NMF factors such that H is close a kinetics matrix.
    [W,H] = scale_factors(W,H);
    
    
    % Plot computed component spectra and relative concentrations.
    if PLOT
        current_folder = sprintf('results\\real_experiment\\r=%i',r);
        if exist(current_folder,'dir') ~= 7
            mkdir(current_folder);
        end
    
        for s = 1:r
            plot(f,W(:,s),'LineWidth',0.6)
            hold on
            yl = ylim;
            for j = 1:r
                line([f(char_idx(j)) f(char_idx(j))], [yl(1),yl(2)],'LineWidth',0.6, 'LineStyle', '-', 'Color', [0, 0, 0, 0.4]);
            end
            hold off
            axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
            fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
            set(gca, axopts{:});
            set(gcf, fopts{:});
            xlabel('Wavenumber (cm^{-1})')%,'FontSize', 18);
            ylabel('Intensity')%,'FontSize', 18);
            xlim([min(f) max(f)])
            print(gcf,strcat(current_folder,'\w',int2str(s),'.png'),'-dpng','-r300'); 
            pause
        end


        
        plot(tspan,H','LineWidth',0.6)
        axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
        fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
        set(gca, axopts{:});
        set(gcf, fopts{:});
        xlabel('Time (s)')%,'FontSize', 18);
        ylabel('Relative concentration')%,'FontSize', 18);
        xlim([min(tspan) max(tspan)])
        ylim([0 1])
        legend_names = {'Species A','Species B','Species C','Species D','Species E','Species F','Species G','Species H'};
        legend_names = {legend_names{1:r}};
        legend(legend_names)%,'Location','northwest')
        print(gcf,strcat(current_folder,'\h.png'),'-dpng','-r300');
        pause
    end

    % Compute reaction coefficients K.
    K = compute_K(H,tspan);

    % Plot comparison of the kinetics based on the computed reaction
    % coefficients and the true kinetics.
    if PLOT
        H_K = get_H(K,tspan,H(:,1)/norm(H(:,1),1));
        plot(tspan,H_K','LineWidth',0.6)
        axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
        fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
        set(gca, axopts{:});
        set(gcf, fopts{:});
        xlabel('Time (s)')%,'FontSize', 18);
        ylabel('Relative concentration')%,'FontSize', 18);
        xlim([min(tspan) max(tspan)])
        ylim([0 1])
        legend(legend_names)%,'Location','northwest')
        print(gcf,strcat(current_folder,'\k.png'),'-dpng','-r300'); 
        pause
    end
    
    % Store results.
    error(counter) = norm(A_without_buffer - W*H,'fro')/norm(A_without_buffer,'fro');
    fprintf('Relative approximation error of data:  %.3e\n', error(counter));
    W_comp{counter} = W;
    H_comp{counter} = H;
    K_comp{counter} = K;
    counter = counter + 1;
end

end
