

%{
f = figure('Visible','off');
%subplot(2,1,1);
plot([0,0,0,6,6,25,25], [0.5045, 0.5226, 0.5350, 0.5422, 0.5647, 0.6681, 0.6718],'x','LineWidth',2);
legend({'$\langle$ dihedral order $\neq 4\rangle $'},'Interpreter','latex')
xlabel('number of fibers','Interpreter','latex');
exportgraphics(f,'dihedral_angle_order_mean.png');
%}

[f, p]     = uigetfile('*.mat');
eval(sprintf('cd %s',p));
true_name  = regexprep(f,'.mat','','ignorecase');
mkdir(true_name);
filename   = [p f];
fileID     = fopen(filename,'r');
load(f);
cleaned_selected_cell = {};
fig_first_clean = figure('Visible','off');
hold on
mean_z_fiber = [];
std_z_fiber  = [];
for i = 1:numel(save_selected_cell)
    cleaned_x = [];
    cleaned_y = [];
    cleaned_z = [];
    for j = 1:2000
        cache_cell_y = [];
        cache_cell_z = [];
        for k = 1:numel(save_selected_cell{i,1})/3
            if save_selected_cell{i,1}(k,1) == j
                cache_cell_y(end+1) = save_selected_cell{i,1}(k,2);
                cache_cell_z(end+1) = save_selected_cell{i,1}(k,3);
            end
        end
        if numel(cache_cell_y) ~=0
            cleaned_x(end+1) = j;
            cleaned_y(end+1) = mean(cache_cell_y);
            cleaned_z(end+1) = mean(cache_cell_z);
        end
    end
    cleaned_x = cleaned_x.';
    cleaned_y = cleaned_y.';
    cleaned_z = cleaned_z.';
    cleaned_selected_cell{end+1} = [cleaned_x,cleaned_y,cleaned_z];
    
    plot3(cleaned_x,cleaned_y,cleaned_z)
    if numel(cleaned_z) ~=0
        mean_z_fiber(end+1) = mean(cleaned_z);
        std_z_fiber(end+1) = std(cleaned_z);
        %plot3([cleaned_x(1),cleaned_x(numel(cleaned_x)-1)],[cleaned_y(1),cleaned_y(numel(cleaned_y)-1)],[mean(cleaned_z),mean(cleaned_z)],'LineWidth',2)
    end
end
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
hold off
view(-45,45)
exportgraphics(fig_first_clean,'3D_fibers_cleaned.png');
view(0,0)
exportgraphics(fig_first_clean,'0_0_3D_fibers_cleaned.png'); 
view(90,0)
exportgraphics(fig_first_clean,'90_0_3D_fibers_cleaned.png'); 
view(0,90)
exportgraphics(fig_first_clean,'0_90_3D_fibers_cleaned.png'); 



%%

figure_z_position_figure = figure('Visible','off');
hold on
errorbar(1:numel(mean_z_fiber),mean_z_fiber,std_z_fiber,'x','LineWidth',2);
for i = 0:4
    if i == 4
        plot([25,29],[mean(mean_z_fiber( 25:29 )), mean(mean_z_fiber( 25:29 ))],'LineWidth',2);
    else
        plot([1+6*i,6+6*i],[mean(mean_z_fiber( (1+6*i):(6+6*i) )), mean(mean_z_fiber( (1+6*i):(6+6*i) ))],'LineWidth',2);
        
    end
    
end
legend({'$z$ position',[num2str(round(mean(mean_z_fiber( 1:6 )))),'$\pm$',num2str(round(mean(std_z_fiber)))],[num2str(round(mean(mean_z_fiber( 7:12 )))),'$\pm$',num2str(round(mean(std_z_fiber)))],[num2str(round(mean(mean_z_fiber( 13:18 )))),'$\pm$',num2str(round(mean(std_z_fiber)))],[num2str(round(mean(mean_z_fiber( 19:24 )))),'$\pm$',num2str(round(mean(std_z_fiber)))],[num2str(round(mean(mean_z_fiber( 25:29 )))),'$\pm$',num2str(round(mean(std_z_fiber)))]},'Location','southwest','Interpreter','Latex');

hold off
xlabel('Associated fiber labelling number','Interpreter','latex');
ylabel('$z$ position (px)','Interpreter','latex')
exportgraphics(figure_z_position_figure,'z_position_fibers.png'); 

%%
%{
load('normalised_quantity_density.mat');
load('zmv.mat');
histogram_z_position_fibers = figure;%('Visible','off');
z_hist = [];
hold on
for i = 1:numel(save_selected_cell)
    cleaned_x = [];
    cleaned_y = [];
    %cleaned_z = [];
    for j = 1:2000
        cache_cell_y = [];
        cache_cell_z = [];
        for k = 1:numel(save_selected_cell{i,1})/3
            if save_selected_cell{i,1}(k,1) == j
                cache_cell_y(end+1) = save_selected_cell{i,1}(k,2);
                cache_cell_z(end+1) = save_selected_cell{i,1}(k,3);
            end
        end
        if numel(cache_cell_y) ~=0
            cleaned_x(end+1) = j;
            cleaned_y(end+1) = mean(cache_cell_y);
            z_hist(end+1) = mean(cache_cell_z);
        end
    end
    cleaned_x = cleaned_x.';
    cleaned_y = cleaned_y.';
    z_hist = z_hist.';

    
end
yyaxis right
ylabel('Count','Interpreter','Latex')
histogram(z_hist,26,'FaceAlpha',0.2);
yyaxis left
ylabel('Number of vertex of order $\neq 4$ divided by total number of vertex','Interpreter','Latex')
errorbar([197.7812,279.5963,369.6744,465.5747,576.5098], [0,0,0,0,0],[18.4692,18.4692,18.4692,18.4692,18.4692],'horizontal','x','LineWidth',2);
plot(zmv,normalised_quantity_density,'x-','LineWidth',2); 
hold off
xlabel('$z$ position [px]','Interpreter','Latex')
exportgraphics(histogram_z_position_fibers,'histogram_position_fibers.png');

%}

%%

load('fibers_cleaned_isolated.mat');

fig_selected = figure('Visible','off');
cleaned_selected_cell = cleaned_selected_cell.';
mean_z_fiber_cleaned = [];
std_z_fiber_cleaned = [];
hold on
for i = 1:numel(cleaned_selected_cell)
    plot_cache_x = [];
    plot_cache_y = [];
    plot_cache_z = [];
    for j = 1:numel(cleaned_selected_cell{i,1})/3
        plot_cache_x(end+1) = cleaned_selected_cell{i,1}(j,1);
        plot_cache_y(end+1) = cleaned_selected_cell{i,1}(j,2);
        plot_cache_z(end+1) = cleaned_selected_cell{i,1}(j,3);
    end
    plot3(plot_cache_x,plot_cache_y,plot_cache_z,'x')
    mean_z_fiber_cleaned(end+1) = mean(plot_cache_z);
    std_z_fiber_cleaned(end+1)  = std(plot_cache_z);
end
std_z_fiber_cleaned(numel(std_z_fiber_cleaned))=std_z_fiber_cleaned(numel(std_z_fiber_cleaned)-1);
hold off
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
%xlim([0, max(save_x)]);
%ylim([0, max(save_y)]);
%zlim([0, max(save_z)]);
view(-45,45)
exportgraphics(fig_selected,'3D_fibers_c.png');
view(0,0)
exportgraphics(fig_selected,'0_0_3D_fibers_c.png'); 
view(90,0)
exportgraphics(fig_selected,'90_0_3D_fibers_c.png'); 
view(0,90)
exportgraphics(fig_selected,'0_90_3D_fibers_c.png'); 

figure_z_cleaned_position_figure = figure('Visible','off');
hold on
errorbar(1:numel(mean_z_fiber_cleaned),mean_z_fiber_cleaned,std_z_fiber_cleaned,'x','LineWidth',2);
for i = 0:4
    if i == 4
        plot([25,29],[mean(mean_z_fiber_cleaned( 25:29 )), mean(mean_z_fiber_cleaned( 25:29 ))],'LineWidth',2);
    else
        plot([1+6*i,6+6*i],[mean(mean_z_fiber_cleaned( (1+6*i):(6+6*i) )), mean(mean_z_fiber_cleaned( (1+6*i):(6+6*i) ))],'LineWidth',2);
        
    end
    
end
legend({'$z$ position',[num2str(round(mean(mean_z_fiber_cleaned( 1:6 )))),'$\pm$',num2str(round(mean(std_z_fiber_cleaned)))],[num2str(round(mean(mean_z_fiber_cleaned( 7:12 )))),'$\pm$',num2str(round(mean(std_z_fiber_cleaned)))],[num2str(round(mean(mean_z_fiber_cleaned( 13:18 )))),'$\pm$',num2str(round(mean(std_z_fiber_cleaned)))],[num2str(round(mean(mean_z_fiber_cleaned( 19:24 )))),'$\pm$',num2str(round(mean(std_z_fiber_cleaned)))],[num2str(round(mean(mean_z_fiber_cleaned( 25:29 )))),'$\pm$',num2str(round(mean(std_z_fiber_cleaned)))]},'Location','southwest','Interpreter','Latex');

hold off
xlabel('Associated fiber labelling number','Interpreter','latex');
ylabel('$z$ position (px)','Interpreter','latex')
exportgraphics(figure_z_cleaned_position_figure,'z_cleaned_position_fibers.png'); 


%%
%{
cleaned_selected_cell2 = {};
fig_first_clean = figure('Visible','off');
hold on
mean_z_skel_fiber = [];
std_z_skel_fiber  = [];
for i = 1:numel(cleaned_selected_cell)
    cleaned_x = [];
    cleaned_y = [];
    cleaned_z = [];
    for j = 1:2000
        cache_cell_y = [];
        cache_cell_z = [];
        for k = 1:numel(cleaned_selected_cell{i,1})/3
            if cleaned_selected_cell{i,1}(k,1) == j
                cache_cell_y(end+1) = cleaned_selected_cell{i,1}(k,2);
                cache_cell_z(end+1) = cleaned_selected_cell{i,1}(k,3);
            end
        end
        if numel(cache_cell_y) ~=0
            cleaned_x(end+1) = j;
            cleaned_y(end+1) = mean(cache_cell_y);
            cleaned_z(end+1) = mean(cache_cell_z);
        end
    end
    cleaned_x = cleaned_x.';
    cleaned_y = cleaned_y.';
    cleaned_z = cleaned_z.';
    cleaned_selected_cell2{end+1} = [cleaned_x,cleaned_y,cleaned_z];
    
    plot3(cleaned_x,cleaned_y,cleaned_z)
    if numel(cleaned_z) ~=0
        mean_z_skel_fiber(end+1) = mean(cleaned_z);
        std_z_skel_fiber(end+1) = std(cleaned_z);
    end
end
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
hold off
view(-45,45)
exportgraphics(fig_first_clean,'3D_fibers_cleaned.png');
view(0,0)
exportgraphics(fig_first_clean,'0_0_3D_fibers_cleaned.png'); 
view(90,0)
exportgraphics(fig_first_clean,'90_0_3D_fibers_cleaned.png'); 
view(0,90)
exportgraphics(fig_first_clean,'0_90_3D_fibers_cleaned.png'); 

%}
%%

read_scale('SLiceRight_Fiber1363.tif');