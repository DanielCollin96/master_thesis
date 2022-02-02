% This script shows how to create some of the figures used in the master thesis
% "Nonnegative matrix factorization and its application to Raman spectral
% data analysis" by Daniel Collin, Technische Universität Berlin. It can be
% modified arbitrarily to create the desired figures.

if exist('figures','dir') ~= 7
    mkdir('figures');
end

% Generate data.
[A, W, H, K, tspan, freq] = generate_real_substance_data();
% A. trans_dichlorethen
% B. toluene
% C. diethylether
% D. cis-dichlorethen
% E. acetonitrile

%% Plot data set.
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
print(gcf,'figures\raman_time_resolved.png','-dpng','-r300'); 



%% Plot Lorentz functions.
lorentz = eval_lorentz_sum([0, 1, 1],linspace(-10,10,1000));

plot(linspace(-10,10,1000),lorentz,'LineWidth',0.6)
axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
set(gca, axopts{:});
set(gcf, fopts{:});
xlabel('x')%,'FontSize', 18);
ylabel('L(x)')%,'FontSize', 18);
xlim([-10 10])
ylim([0 1.1])
print(gcf,'figures\lorentzian1.png','-dpng','-r300'); 


lorentz = eval_lorentz_sum([1950, 5, 30000],freq);

plot(freq,lorentz,'LineWidth',0.6)
axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
set(gca, axopts{:});
set(gcf, fopts{:});
xlabel('x')%,'FontSize', 18);
ylabel('L(x)')%,'FontSize', 18);
xlim([min(freq) max(freq)])
ylim([0 35000])
print(gcf,'figures\lorentzian2.png','-dpng','-r300');



%% Plot comparison of Lorentz function fit and real data.
param = greedy_lorentz(freq,W(:,2),0);
lorentz = eval_lorentz_sum(param,freq);

plot(freq,lorentz,'LineWidth',0.6)
hold on
plot(freq,W(:,2),'k--','LineWidth',0.6)
hold off
axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
set(gca, axopts{:});
set(gcf, fopts{:});
xlabel('x')%,'FontSize', 18);
ylabel('w_s(x)')%,'FontSize', 18);
xlim([min(freq) max(freq)])
%ylim([0 35000])
legend('Fitted Lorentzians','Original data')
print(gcf,'figures\lorentz_fit1.png','-dpng','-r300');


param = greedy_lorentz(freq,W(:,5),0);
lorentz = eval_lorentz_sum(param,freq);

plot(freq,lorentz,'LineWidth',0.6)
hold on
plot(freq,W(:,5),'k--','LineWidth',0.6)
hold off
axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
set(gca, axopts{:});
set(gcf, fopts{:});
xlabel('x')%,'FontSize', 18);
ylabel('w_s(x)')%,'FontSize', 18);
xlim([min(freq) max(freq)])
%ylim([0 35000])
legend({'Fitted Lorentzians','Original data'},'Location','northwest')
print(gcf,'figures\lorentz_fit2.png','-dpng','-r300');



%% Plot component spectra and relative concentrations.
rng(10)
n_species = 4;
freq = linspace(500,3500,2501);

p = generate_rand_params(freq,n_species);

W = zeros(length(freq),n_species);
for j = 1:n_species
    W(:,j) = eval_lorentz_sum(p{j},freq);
end

plot(freq,W,'LineWidth',0.6)
axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
set(gca, axopts{:});
set(gcf, fopts{:});
xlabel('Wavenumber (cm^{-1})')%,'FontSize', 18);
ylabel('w_s')%,'FontSize', 18);
xlim([min(freq) max(freq)])
%ylim([0 max(max(W))*1.01])
legend({'s = 1','s = 2','s = 3','s = 4'},'Location','northwest')
print(gcf,'figures\w_example.png','-dpng','-r300'); 

plot(tspan,H','LineWidth',0.6)
axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
set(gca, axopts{:});
set(gcf, fopts{:});
xlabel('t')%,'FontSize', 18);
ylabel('h_s(t)')%,'FontSize', 18);
xlim([min(tspan) max(tspan)])
%ylim([0 max(max(W))*1.01])
legend({'Species A','Species B','Species C','Species D','Species E'},'Location','northwest')
print(gcf,'figures\h_example.png','-dpng','-r300'); 


%% Plot example 3.6 from master thesis.
rng(10)
n_species = 3;
freq = linspace(500,3500,2501);
tspan = linspace(0,20,151);

p = generate_rand_params(freq,n_species);

W = zeros(length(freq),n_species);
for j = 1:n_species
    W(:,j) = eval_lorentz_sum(p{j},freq);
end

plot(freq,W,'LineWidth',0.6)
axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
set(gca, axopts{:});
set(gcf, fopts{:});
xlabel('Wavenumber (cm^{-1})')%,'FontSize', 18);
ylabel('Intensity')%,'FontSize', 18);
xlim([min(freq) max(freq)])
%ylim([0 max(max(W))*1.01])
legend({'Species A','Species B','Species C'},'Location','northwest')
print(gcf,'figures\raman_example_w.png','-dpng','-r300'); 

K = [-0.5   0 0;
      0.4 -0.25 0;
      0.1  0.25 0];
h0 = [1 0 0]';
H = get_H(K,tspan,h0);

plot(tspan,H','LineWidth',0.6)
axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
set(gca, axopts{:});
set(gcf, fopts{:});
xlabel('Time (s)')%,'FontSize', 18);
ylabel('Relative concentration')%,'FontSize', 18);
xlim([min(tspan) max(tspan)])
%ylim([0 max(max(W))*1.01])
legend({'Species A','Species B','Species C'},'Location','northwest')
print(gcf,'figures\raman_example_h.png','-dpng','-r300'); 


A = W*H;
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
print(gcf,'figures\raman_example_a.png','-dpng','-r300'); 

raman_analysis(3,A,W,H,K,tspan,freq,0.01,2);


%% Plot examples of component spectra with different levels of interference and their characteristic frequencies.
rng(3)
n_species = 5;
freq = linspace(500,3500,2501);

p = generate_rand_params(freq,n_species);

W = zeros(length(freq),n_species);
for j = 1:n_species
    W(:,j) = eval_lorentz_sum(p{j},freq);
end

% Compute characteristic frequencies.
W(~any(W,2),:)= [];
[m,r] = size(W);
W_norm = zeros(m,r);
for i = 1:m
    W_norm(i,:) = W(i,:)/norm(W(i,:),1);
end
val = zeros(r,1);
idx = zeros(r,1);
for i = 1:r
    [val(i),idx(i)] = max(W_norm(:,i));
end

plot(freq,W,'LineWidth',0.6)
hold on
yl = ylim;
for i = 1:r
    line([freq(idx(i)) freq(idx(i))], [yl(1),yl(2)], 'LineWidth',0.6, 'LineStyle', '-', 'Color', [0, 0, 0, 0.4]);
end
hold off
axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
set(gca, axopts{:});
set(gcf, fopts{:});
xlabel('Wavenumber (cm^{-1})')%,'FontSize', 18);
ylabel('Intensity')%,'FontSize', 18);
xlim([min(freq) max(freq)])
%ylim([0 max(max(W))*1.01])
legend({'Species A','Species B','Species C','Species D','Species E'},'Location','northwest')
print(gcf,'figures\interference_example1.png','-dpng','-r300'); 



rng(8)
n_species = 3;
freq = linspace(500,3500,2501);
n_meet = 1;

p = generate_rand_params(freq,n_species);
p = move_peaks(p,freq,n_meet,0.8);

W = zeros(length(freq),n_species);
for j = 1:n_species
    W(:,j) = eval_lorentz_sum(p{j},freq);
end

% Compute characteristic frequencies.
W(~any(W,2),:)= [];
[m,r] = size(W);

W_norm = zeros(m,r);
for i = 1:m
    W_norm(i,:) = W(i,:)/norm(W(i,:),1);
end

val = zeros(r,1);
idx = zeros(r,1);
for i = 1:r
    [val(i),idx(i)] = max(W_norm(:,i));
end

plot(freq,W,'LineWidth',0.6)
hold on
yl = ylim;
for i = 1:r
    line([freq(idx(i)) freq(idx(i))], [yl(1),yl(2)], 'LineWidth',0.6, 'LineStyle', '-', 'Color', [0, 0, 0, 0.4]);
end
hold off
axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1.2 1 1]};%,'FontSize', 18};
fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
set(gca, axopts{:});
set(gcf, fopts{:});
xlabel('Wavenumber (cm^{-1})')%,'FontSize', 18);
ylabel('Intensity')%,'FontSize', 18);
xlim([min(freq) max(freq)])
%ylim([0 max(max(W))*1.01])
legend({'Species A','Species B','Species C'},'Location','northwest')
print(gcf,'figures\interference_example2.png','-dpng','-r300');


