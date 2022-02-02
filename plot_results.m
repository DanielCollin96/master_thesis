function plot_results(exp,data_file)
% This function plots selected results of the numerical experiment on the
% effects of increasing interference.

% Load results.
data = load(strcat('results\',data_file));
result = data.result;

if exp == 1 || exp == 2
    noise = result.noise;
else
    noise = result{1}.noise;
end



% Create folder for figures.
filename = strcat('results\exp_',num2str(exp),'_noise_',num2str(noise));

if exist(filename,'dir') ~= 7
    mkdir(filename);
end


if exp == 1 || exp == 2
    
    average = result.average;
    average_interference = result.average_interference;
    
    % Plot average relative approximation error of A.
    plot(average_interference,average.SPA(1,:),'x-','LineWidth',0.6);
    hold on
    plot(average_interference,average.SNPA(1,:),'o-','LineWidth',0.6);
    plot(average_interference,average.XRAY(1,:),'+-','LineWidth',0.6);
    plot(average_interference,average.ALS(1,:),'s-','LineWidth',0.6);
    axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1 0.6 1]};%,'FontSize', 18};
    fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
    set(gca, axopts{:});
    set(gcf, fopts{:});
    title('Average relative error E_A');
    xlabel('I(W)');
    ylabel('Error');
    %xlim([0 1])
    %ylim([0 1])
    legend({'SPA', 'SNPA', 'XRAY','ALS'},'Location','northwest'); 
    print(gcf,strcat(filename,'\a.png'),'-dpng','-r300'); 
    hold off

    % Plot average relative approximation error of W.
    plot(average_interference,average.SPA(2,:),'x-','LineWidth',0.6);
    hold on
    plot(average_interference,average.SNPA(2,:),'o-','LineWidth',0.6);
    plot(average_interference,average.XRAY(2,:),'+-','LineWidth',0.6);
    plot(average_interference,average.ALS(2,:),'s-','LineWidth',0.6);
    axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1 0.6 1]};%,'FontSize', 18};
    fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
    set(gca, axopts{:});
    set(gcf, fopts{:});
    title('Average relative error E_W');
    xlabel('I(W)');
    ylabel('Error');
    %xlim([0 1])
    %ylim([0 3])
    legend({'SPA', 'SNPA', 'XRAY','ALS'},'Location','northwest'); 
    print(gcf,strcat(filename,'\w.png'),'-dpng','-r300'); 
    hold off

    % Plot average relative approximation error of H.
    plot(average_interference,average.SPA(3,:),'x-','LineWidth',0.6);
    hold on
    plot(average_interference,average.SNPA(3,:),'o-','LineWidth',0.6);
    plot(average_interference,average.XRAY(3,:),'+-','LineWidth',0.6);
    plot(average_interference,average.ALS(3,:),'s-','LineWidth',0.6);
    axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1 0.6 1]};%,'FontSize', 18};
    fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
    set(gca, axopts{:});
    set(gcf, fopts{:});
    title('Average relative error E_H');
    xlabel('I(W)');
    ylabel('Error');
    %xlim([0 1])
    %ylim([0 1])
    legend({'SPA', 'SNPA', 'XRAY','ALS'},'Location','northwest'); 
    print(gcf,strcat(filename,'\h.png'),'-dpng','-r300'); 
    hold off

    % Plot average relative approximation error of K.
    plot(average_interference,average.SPA(4,:),'x-','LineWidth',0.6);
    hold on
    plot(average_interference,average.SNPA(4,:),'o-','LineWidth',0.6);
    plot(average_interference,average.XRAY(4,:),'+-','LineWidth',0.6);
    plot(average_interference,average.ALS(4,:),'s-','LineWidth',0.6);
    axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1 0.6 1]};%,'FontSize', 18};
    fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
    set(gca, axopts{:});
    set(gcf, fopts{:});
    title('Average relative error E_K');
    xlabel('I(W)');
    ylabel('Error');
    %xlim([0 1])
    %ylim([0 1])
    legend({'SPA', 'SNPA', 'XRAY','ALS'},'Location','northwest'); 
    print(gcf,strcat(filename,'\k.png'),'-dpng','-r300'); 
    hold off
    
else
    
    % For experiment 3:
    
    x = zeros(length(result{1}.average_interference),length(result));
    y = zeros(length(result{1}.average_interference),length(result));
    legend_names = cell(length(result),1);
    for i = 1:length(result)
        x(:,i) = result{i}.average_interference';
        y(:,i) = result{i}.average.SNPA(1,:)';
        legend_names{i} = strcat('r = ',num2str(result{i}.r));
    end

    % Plot average relative approximation error of A.
    plot(x,y,'x-','LineWidth',0.6);
    axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1 0.6 1]};%,'FontSize', 18};
    fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
    set(gca, axopts{:});
    set(gcf, fopts{:});
    title('Average relative error E_A');
    xlabel('I(W)');
    ylabel('Error');
    %xlim([0 1])
    %ylim([0 1])
    legend(legend_names,'Location','northwest'); 
    print(gcf,strcat(filename,'\a.png'),'-dpng','-r300'); 

    % Plot average relative approximation error of W.
    for i = 1:length(result)
        y(:,i) = result{i}.average.SNPA(2,:)';
    end
    plot(x,y,'x-','LineWidth',0.6);
    axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1 0.6 1]};%,'FontSize', 18};
    fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
    set(gca, axopts{:});
    set(gcf, fopts{:});
    title('Average relative error E_W');
    xlabel('I(W)');
    ylabel('Error');
    %xlim([0 1])
    %ylim([0 1])
    legend(legend_names,'Location','northwest'); 
    print(gcf,strcat(filename,'\w.png'),'-dpng','-r300'); 


    % Plot average relative approximation error of H.
    for i = 1:length(result)
        y(:,i) = result{i}.average.SNPA(3,:)';
    end
    plot(x,y,'x-','LineWidth',0.6);
    axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1 0.6 1]};%,'FontSize', 18};
    fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
    set(gca, axopts{:});
    set(gcf, fopts{:});
    title('Average relative error E_H');
    xlabel('I(W)');
    ylabel('Error');
    %xlim([0 1])
    %ylim([0 1])
    legend(legend_names,'Location','northwest'); 
    print(gcf,strcat(filename,'\h.png'),'-dpng','-r300'); 

    % Plot average relative approximation error of K.
    for i = 1:length(result)
        y(:,i) = result{i}.average.SNPA(4,:)';
    end
    plot(x,y,'x-','LineWidth',0.6);
    axopts = {'XMinorTick', 'on', 'PlotBoxAspectRatio', [1 0.6 1]};%,'FontSize', 18};
    fopts = {'Position', [0 0 560 420], 'PaperPositionMode', 'auto'};
    set(gca, axopts{:});
    set(gcf, fopts{:});
    title('Average relative error E_K');
    xlabel('I(W)');
    ylabel('Error');
    %xlim([0 1])
    %ylim([0 1])
    legend(legend_names,'Location','northwest'); 
    print(gcf,strcat(filename,'\k.png'),'-dpng','-r300'); 

    
end

end