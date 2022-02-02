function plot_mask(A, mask, freq, heur)
% This function visualizes a selection of freqencies, specified as mask.
% Depending on which heuristic is selected, the properties of the plot
% differ slightly. 
% heur = 1: beta heuristic
% heur = 2: gamma heuristic
% heur = 3: both heuristics combined

m = size(A,1);

% Get maximum and minumum intensity of each frequency.
maxA = max(A, [], 2);
minA = min(A, [], 2);

% Get first and last spectrum of data.
first = A(:,1);
last = A(:,end);

% Plot characterisics of data set.
plot(freq, [first, last, minA, maxA],'LineWidth',0.6);
hold on

% Get selected indices.
selection = find(mask == 1);
n_selected = length(selection);

% Get upper and lower bound of vertical direction of plot.
yl = ylim;

% Plot selected frequencies as horizontal lines.
for i = 1:n_selected
    idx = selection(i);
    if heur == 1
        line([freq(idx) freq(idx)], [yl(1),yl(2)], 'LineStyle', '-', 'Color', [0.5, 0.5, 0.5, 0.1]);
    elseif heur == 2
        line([freq(idx) freq(idx)], [yl(1),yl(2)], 'LineStyle', '-', 'Color', [0.5, 0.5, 0.5, 0.3]);
    elseif heur == 3
        line([freq(idx) freq(idx)], [yl(1),yl(2)], 'LineWidth', 0.6, 'LineStyle', '-', 'Color', [0, 0, 0, 0.4]);
    end
end

% Set properties of plot.
axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
set(gca, axopts{:});
set(gcf, fopts{:});
xlabel('Wavenumber (cm^{-1})')%,'FontSize', 18);
ylabel('Intensity')%,'FontSize', 18);
xlim([min(freq) max(freq)])
legend({'First','Last','Min','Max'},'Location','northwest')

% Save figure.
if heur == 1
    print(gcf,'figures\alg_example_beta_heuristic.png','-dpng','-r300');
elseif heur == 2
    print(gcf,'figures\alg_example_gamma_heuristic.png','-dpng','-r300');
elseif heur == 3
    print(gcf,'figures\alg_example_combined_heuristic.png','-dpng','-r300');
end

hold off

fprintf('Number of selected frequencies: %i. Number of removed frequencies: %i.\n',n_selected, m-n_selected);

end