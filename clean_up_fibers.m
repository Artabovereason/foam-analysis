close all
clear
clc

load('fibers_isolated_SlicesRight_Fiber.mat');

cleaned_selected_cell = {};
number_fibers=30;
for n = 1:30
    save_x = [];
    save_y = [];
    save_z = [];
    for k = 1:numel(save_selected_cell{n,1})/3
        save_x(end+1) = save_selected_cell{n,1}(k,1);
        save_y(end+1) = save_selected_cell{n,1}(k,2);
        save_z(end+1) = save_selected_cell{n,1}(k,3);
    end
    figure_to_check = figure;
    plot3(save_x,save_y,save_z,'x');
    view(-45,45)
    %exportgraphics(figure_to_check,'un-cleaned_fiber.png');

    pause(2);
    nb_regions=input('How many region do you want to isolate ?');
    nb_regions = 2*nb_regions;
    figure_to_tag = figure;
    plot3(save_x,save_y,save_z,'x');
    view(0,0)
    pointts_0_0 = ginput3d(nb_regions);
    view(90,0)
    pointts_90_0 = ginput3d(nb_regions);
    view(0,90)
    pointts_0_90 = ginput3d(nb_regions);

    save_selected_x = [];
    save_selected_y = [];
    save_selected_z = [];
    for w = 1:nb_regions/2
        for i=1:numel(save_x)
            if save_x(i) > pointts_0_0(2*w,1) && save_x(i) < pointts_0_0(-1+2*w,1)   
                if save_z(i) > pointts_0_0(2*w,3) && save_z(i) < pointts_0_0(-1+2*w,3)
                    if save_y(i) > pointts_90_0(2*w,2) && save_y(i) < pointts_90_0(-1+2*w,2)   
                        if save_z(i) > pointts_90_0(2*w,3) && save_z(i) < pointts_90_0(-1+2*w,3)
                            if save_x(i) > pointts_0_90(2*w,1) && save_x(i) < pointts_0_90(-1+2*w,1)   
                                if save_y(i) > pointts_0_90(2*w,2) && save_y(i) < pointts_0_90(-1+2*w,2)   
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
    end
    
    save_selected_x = save_selected_x.';
    save_selected_y = save_selected_y.';
    save_selected_z = save_selected_z.';
    
    index_to_delete = [];
    for p = 1:numel(save_x)
        for k = 1:numel(save_selected_x)
            if save_x(p) == save_selected_x(k) && save_y(p) == save_selected_y(k)  && save_z(p) == save_selected_z(k) 
                index_to_delete(end+1) = p;
            end
        end
    end
    cleaned_selected_x = [];
    cleaned_selected_y = [];
    cleaned_selected_z = [];
    for p = 1:numel(save_x)
        cache = 0;
        for  k = 1:numel(index_to_delete)
            if p == index_to_delete(k)
                cache = 1;
            end
        end
        if cache == 0
            cleaned_selected_x(end+1) = save_x(p);
            cleaned_selected_y(end+1) = save_y(p);
            cleaned_selected_z(end+1) = save_z(p);
        end
    end
    cleaned_selected_x = cleaned_selected_x.';
    cleaned_selected_y = cleaned_selected_y.';
    cleaned_selected_z = cleaned_selected_z.';
    
    cleaned_selected_cell{end+1} = [cleaned_selected_x,cleaned_selected_y,cleaned_selected_z];
    
    
    figure_to_sure = figure;
    plot3(cleaned_selected_x,cleaned_selected_y,cleaned_selected_z,'x');
    view(-45,45)
    %exportgraphics(figure_to_sure,'cleaned_fiber.png');
end
%save('fibers_cleaned_isolated.mat','cleaned_selected_cell');

