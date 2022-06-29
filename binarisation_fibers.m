close all
clear
clc
gamma_26   = 0.026;
ext        = '*.tif';
count      = 0;
Nom_essai  = {zeros};
dF_cell    = {zeros};
l_cell     = {zeros};
abs_cell   = {zeros};
ord_cell   = {zeros};
moy_abs    = {zeros};
moy_ord    = {zeros};
err_abs    = {zeros};
err_ord    = {zeros};

%% Treatment
[fi, rep]     = uigetfile('*.tif');
chemin        = fullfile(rep, '*.tif');
list          = dir(chemin); 
save_all_data = {};
    
for w = 170:1480%numel(list)
    count=count+1;
    %% RESIZE TOMOGRAPH IMAGES
    image = fullfile(rep, list(w).name);
    I     = imread(image);       %Original image
    gs    = rgb2gray(I);         %Shading of grey
    I_BW  = imbinarize(gs,0.35); %Binarization
     
    %imshow(I_BW);
    %disp(size(I_BW));
    %[xi,yi] = ginputWhite(2); %Using the command 'ginputWhite' to define a ROI
    %[xi,yi] = ginputWhite(2); %Using the command 'ginputWhite' to define a ROI
    x1 = 0;
    y1 = 0;
    x2 = max(size(I_BW));
    y2 = min(size(I_BW));
    
    for i = 1:min(size(I_BW))
        for j = 1:max(size(I_BW))
            if i < 20 
                I_BW(i,j) = 0;
            end 
            if i > min(size(I_BW))-20 
                I_BW(i,j) = 0;
            end 
            
            if j < 20 
                I_BW(i,j) = 0;
            end 
            if i > max(size(I_BW))-125 
                I_BW(i,j) = 0;
            end 

        end
    end
    cache = {};
    for i = 1:min(size(I_BW))
        for j = 1:max(size(I_BW))
            if I_BW(i,j) ~= 0
                cache{end+1} = [j i];
            end
        end
    end
    save_all_data{end+1} = cache;
end
save_all_data = save_all_data.';
matrix_data   = {};
for i=1:numel(save_all_data)
    for j=1:numel(save_all_data{i,1})
        matrix_data{end+1} = [i,save_all_data{i,1}{1,j}(1),save_all_data{i,1}{1,j}(2)];
    end
end
matrix_data          = matrix_data.';
matrix_data_to_table = zeros(numel(matrix_data),3);
save_x               = zeros(numel(matrix_data),1);
save_y               = zeros(numel(matrix_data),1);
save_z               = zeros(numel(matrix_data),1);
for i = 1:numel(matrix_data)
    matrix_data_to_table(i,1) = matrix_data{i,1}(1);
    matrix_data_to_table(i,2) = matrix_data{i,1}(2);
    matrix_data_to_table(i,3) = matrix_data{i,1}(3);
    save_x(i) = matrix_data{i,1}(1);
    save_y(i) = matrix_data{i,1}(2);
    save_z(i) = matrix_data{i,1}(3);
end
fig3D = figure;
xlim([0, max(save_x)]);
ylim([0, max(save_y)]);
zlim([0, max(save_z)]);
plot3(save_x,save_y,save_z,'x');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
view(-45,45)
exportgraphics(fig3D,'3D_fibers.png');
view(0,0)
exportgraphics(fig3D,'0_0_3D_fibers.png'); 
view(90,0)
exportgraphics(fig3D,'90_0_3D_fibers.png'); 
view(0,90)
exportgraphics(fig3D,'0_90_3D_fibers.png'); 
number_fibers = 30;
save_selected_cell = {};
for n = 1:number_fibers
    plot3(save_x,save_y,save_z,'x');
    view(0,0)
    pointts_0_0  = ginput3d(2);
    view(90,0)
    pointts_90_0 = ginput3d(2);
    view(0,90)
    pointts_0_90 = ginput3d(2);
    save_selected_x = [];
    save_selected_y = [];
    save_selected_z = [];
    
    for i=1:numel(save_x)
        if save_x(i) > pointts_0_0(2,1) && save_x(i) < pointts_0_0(1,1)   
            if save_z(i) > pointts_0_0(2,3) && save_z(i) < pointts_0_0(1,3)
                if save_y(i) > pointts_90_0(2,2) && save_y(i) < pointts_90_0(1,2)   
                    if save_z(i) > pointts_90_0(2,3) && save_z(i) < pointts_90_0(1,3)
                        if save_x(i) > pointts_0_90(2,1) && save_x(i) < pointts_0_90(1,1)   
                            if save_y(i) > pointts_0_90(2,2) && save_y(i) < pointts_0_90(1,2)   
                                save_selected_x(end+1) = save_x(i);
                                save_selected_y(end+1) = save_y(i);
                                save_selected_z(end+1) = save_z(i);
                            end
                        end
                    end
                end   
            end
        end
    end
    save_selected_x = save_selected_x.';
    save_selected_y = save_selected_y.';
    save_selected_z = save_selected_z.';
    save_selected_cell{end+1} = [save_selected_x,save_selected_y,save_selected_z];
    
end
%%
fig_selected = figure('Visible','off');
save_selected_cell = save_selected_cell.';
save('fibers_isolated.mat','save_selected_cell');
hold on
for i = 1:numel(save_selected_cell)
    plot_cache_x = [];
    plot_cache_y = [];
    plot_cache_z = [];
    for j = 1:numel(save_selected_cell{i,1})/3
        plot_cache_x(end+1) = save_selected_cell{i,1}(j,1);
        plot_cache_y(end+1) = save_selected_cell{i,1}(j,2);
        plot_cache_z(end+1) = save_selected_cell{i,1}(j,3);
    end
    plot3(plot_cache_x,plot_cache_y,plot_cache_z,'x')
end
hold off
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$z$','Interpreter','latex');
xlim([0, max(save_x)]);
ylim([0, max(save_y)]);
zlim([0, max(save_z)]);
view(-45,45)
exportgraphics(fig_selected,'3D_fibers_s.png');
view(0,0)
exportgraphics(fig_selected,'0_0_3D_fibers_s.png'); 
view(90,0)
exportgraphics(fig_selected,'90_0_3D_fibers_s.png'); 
view(0,90)
exportgraphics(fig_selected,'0_90_3D_fibers_s.png'); 

