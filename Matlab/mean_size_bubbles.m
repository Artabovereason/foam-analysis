%% Introduction
% uigetfile    : allows us to pick the file we want, one may want to modify
% '*.csv' for other extensions.
%
% delimiter    : what is defining the delimitation in-between cells.
% startRow     : the row at which the data is, meaning we skip the name of the
% cells.
%
% formatSpec   : all of the data we have here is double, and we have 28
% columns, so we have 28 times '%f'.
%
% dataArray    : it is the array in which all the data from all the columns is. 
%
% size_arrays  : define the length of the arrays.
%
% x/y/z_array  : 'Pos. x' , 'Pos. y' and 'Pos. z' columns of the data.
% volume_array : 'Volume' column of the data.
%
% PCA1_x etc   : PCA column associated.
%
% mean_volume     : calculation of the mean volume.
% variance_volume : calculation of the variance of the volume.
% std_volume      : calculation of the standard deviation of the volume.


[f, p]     = uigetfile('*.csv');
eval(sprintf('cd %s',p));
true_name  = regexprep(f,'.csv','','ignorecase');
mkdir(true_name);
filename   = [p f];
fileID     = fopen(filename,'r');
delimiter  = ';';
startRow   = 2;
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
dataArray  = textscan(fileID,formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false,'EndOfLine', '\r\n');
fclose(fileID);
size_arrays  = size(dataArray{1,1});
size_arrays  = size_arrays(1);

x_array      = dataArray{1,1} ;
y_array      = dataArray{1,2} ;
z_array      = dataArray{1,3} ;
volume_array = dataArray{1,4} ;

sphericity   = dataArray{1,7} ;

PCA1_x       = dataArray{1,9} ;
PCA1_y       = dataArray{1,10};
PCA1_z       = dataArray{1,11};
PCA2_x       = dataArray{1,12};
PCA2_y       = dataArray{1,13};
PCA2_z       = dataArray{1,14};
PCA3_x       = dataArray{1,15};
PCA3_y       = dataArray{1,16};
PCA3_z       = dataArray{1,17};

mean_volume     = mean(volume_array);
variance_volume = var(volume_array);
std_volume      = std(volume_array);
disp('Mean volume :')
disp(mean_volume);
disp('Variance volume :')
disp(variance_volume);
disp('Standard deviation volume :');
disp(std_volume);

%% 3D plot of the center of the cells 
plot3(x_array,y_array,z_array,'o');
xlabel('$x$ position [mm]','Interpreter','Latex');
ylabel('$y$ position [mm]','Interpreter','Latex');
zlabel('$z$ position [mm]','Interpreter','Latex');
title('3D plot of cell''s centers','Interpreter','latex')

%% Plot of the volume as a function of the z position

fig0 = figure('visible','off');
scatter(z_array,volume_array,'x','LineWidth',2);
title('Volume of the bubbles as a function of the $z$ position','Interpreter','Latex');
xlabel('$z$ position [mm]','Interpreter','Latex');
ylabel('Volume [mm$^3$]','Interpreter','Latex');
exportgraphics(fig0,strcat(true_name,'/Volume_fct_of_z.png'));

%% histogram of the volumes distribution
fig_hist = figure('visible','off');
histogram(volume_array,15);
ylabel('Counts','Interpreter','latex');
xlabel('Volume [mm$^3$]','Interpreter','latex');
exportgraphics(fig_hist,strcat(true_name,'/histogram_volume.png'));

%% Binning of z axis to calculate intermediate mean values of the volume
% num_slices    : number of even slices of our z axis
% bin_slice_z   : define the z value limit of each bins
% plot_slice_z  : mean value z of each bins
% mean_volume_z : mean volume of each z-bins
% std_volume_z  : std of volume of each z-bins
% 
num_slices    = 21;
bin_slice_z   = linspace(min(z_array),max(z_array),num_slices);
plot_slice_z  = [];
for k = 1:(numel(bin_slice_z)-1)
    plot_slice_z(end+1) = (bin_slice_z(k)+bin_slice_z(k+1))/2;
end
mean_volume_z = [];
std_volume_z  = [];
num_bin       = zeros(num_slices-1,1);
for j = 1:(num_slices-1)
    cache_array  = [];
    for i = 1:size_arrays
        if z_array(i) > bin_slice_z(j) && z_array(i) < bin_slice_z(j+1)
            cache_array(end+1) = volume_array(i);
            num_bin(j) = num_bin(j)+1;
        end
    end
    mean_volume_z(end+1)= mean(cache_array);
    std_volume_z(end+1) = std(cache_array);
end

num_bin = num_bin./size_arrays;
mean_volume_z = mean_volume_z.'; % .' is doing the transposition of the matrix same as transpose(A)
std_volume_z  = std_volume_z.';
mean_mean_volume_z = 0;
for i=1:numel(mean_volume_z)
    mean_mean_volume_z = mean_mean_volume_z+num_bin(i)*mean_volume_z(i);
end
%mean_mean_volume_z= mean_mean_volume_z/numel(mean_volume_z);


%% Plot

fig1 = figure('visible','off');
error_bar_mean = [] ;
for l = 1:numel(plot_slice_z)
    error_bar_mean(end+1) = plot_slice_z(1)-bin_slice_z(1);
end

errorbar(plot_slice_z,mean_volume_z,error_bar_mean,'horizontal','.','LineWidth',2)

hold on 
plot([min(bin_slice_z),max(bin_slice_z)],[mean(volume_array) , mean(volume_array)]  , 'LineWidth',2)
plot([min(bin_slice_z),max(bin_slice_z)],[mean(mean_volume_z), mean(mean_volume_z)] , 'LineWidth',2)
legend({'$\langle V\rangle(z)$','$\langle \langle V\rangle(z)\rangle$ ','$\langle V\rangle$'},'Location','southwest','Interpreter','Latex');
hold off
xlim([min(bin_slice_z),max(bin_slice_z)])
xlabel('$z$ position [mm]','Interpreter','Latex');
exportgraphics(fig1,strcat(true_name,'/mean_value_volume.png'));

%% Plot

fig2 = figure('visible','off');
%plot(plot_slice_z,std_volume_z,'x','LineWidth',2)
error_bar_std = [] ;
for l = 1:numel(plot_slice_z)
    error_bar_std(end+1) = plot_slice_z(1)-bin_slice_z(1);
end

errorbar(plot_slice_z,std_volume_z,error_bar_std,'horizontal','.','LineWidth',2)
hold on 
plot([min(bin_slice_z),max(bin_slice_z)],[std_volume, std_volume], 'LineWidth',2)
legend({'std$(z)$','std'},'Location','southwest','Interpreter','Latex');
title('Standard Deviation as a function of the $z$ position','Interpreter','Latex')
hold off
xlim([min(bin_slice_z),max(bin_slice_z)])
xlabel('$z$ position [mm]','Interpreter','Latex');
exportgraphics(fig2,strcat(true_name,'/std_z.png'));


%% Plot

fig4 = figure('visible','off');
subplot(2,2,1);
scatter(z_array,volume_array,'x','LineWidth',1);
title('Volume as a function of the $z$ position','Interpreter','Latex');
xlabel('$z$ position [mm]','Interpreter','Latex');
ylabel('Volume [mm$^3$]','Interpreter','Latex');
xlim([min(bin_slice_z),max(bin_slice_z)])

subplot(2,2,2);
hold on
plot([min(bin_slice_z),max(bin_slice_z)],[mean(volume_array) , mean(volume_array)]  , 'LineWidth',1)
%plot([min(bin_slice_z),max(bin_slice_z)],[mean(mean_volume_z), mean(mean_volume_z)] , 'LineWidth',1)
plot([min(bin_slice_z),max(bin_slice_z)],[mean_mean_volume_z , mean_mean_volume_z] ,'--', 'LineWidth',1)

errorbar(plot_slice_z,mean_volume_z,error_bar_mean,'horizontal','.','LineWidth',1,'Color','black')
legend({'$\langle \langle V\rangle(z)\rangle$ ','$\langle V\rangle$','$\langle V\rangle(z)$'},'Location','southwest','Interpreter','Latex');
xlim([min(bin_slice_z),max(bin_slice_z)])
xlabel('$z$ position [mm]','Interpreter','Latex');
hold off

subplot(2,2,3)
hold on
plot([min(bin_slice_z),max(bin_slice_z)],[std_volume, std_volume], 'LineWidth',1)
errorbar(plot_slice_z,std_volume_z,error_bar_std,'horizontal','.','LineWidth',1,'Color','black')
legend({'$\langle$std$(z)\rangle$','Standard Deviation of $z$'},'Location','southwest','Interpreter','Latex');
title('Standard Deviation as a function of the $z$ position','Interpreter','Latex')
xlim([min(bin_slice_z),max(bin_slice_z)])
xlabel('$z$ position [mm]','Interpreter','Latex');

sgtitle(string(f),'Interpreter', 'none') 
hold off

exportgraphics(fig4,strcat(true_name,'/all_statistics.png'));

%% Plot

fig5 = figure('visible','off');
subplot(3,1,1)
scatter(PCA1_x,PCA2_x,'x');
xlabel('PCA1 on $x$ axis','Interpreter','Latex')
ylabel('PCA2 on $x$ axis','Interpreter','Latex')

subplot(3,1,2)
scatter(PCA1_y,PCA2_y,'x');
xlabel('PCA1 on $y$ axis','Interpreter','Latex')
ylabel('PCA2 on $y$ axis','Interpreter','Latex')

subplot(3,1,3)
scatter(PCA1_z,PCA2_z,'x');
xlabel('PCA1 on $z$ axis','Interpreter','Latex')
ylabel('PCA2 on $z$ axis','Interpreter','Latex')

sgtitle(string(f),'Interpreter', 'none') 
exportgraphics(fig5,strcat(true_name,'/PCA_tests.png'));

%% Plot
%{
fig6 = figure;
plot3(PCA1_x,PCA2_x,PCA3_x,'.');
xlabel('PCA1 on $x$ axis','Interpreter','Latex')
ylabel('PCA2 on $x$ axis','Interpreter','Latex')
zlabel('PCA3 on $x$ axis','Interpreter','Latex')

fig7 = figure;
plot3(PCA1_y,PCA2_y,PCA3_y,'.');
xlabel('PCA1 on $y$ axis','Interpreter','Latex')
ylabel('PCA2 on $y$ axis','Interpreter','Latex')
zlabel('PCA3 on $y$ axis','Interpreter','Latex')


fig8 = figure;
plot3(PCA1_z,PCA2_z,PCA3_z,'.');
xlabel('PCA1 on $z$ axis','Interpreter','Latex')
ylabel('PCA2 on $z$ axis','Interpreter','Latex')
zlabel('PCA3 on $z$ axis','Interpreter','Latex')

%}
 
%% pair correlation function calculation 
% total_sample_v : total volume of the sample
% density_sample : density of the sample

total_sample_v = (max(x_array)-min(x_array)).*(max(y_array)-min(y_array)).*(max(z_array)-min(z_array));
density_sample = size_arrays/total_sample_v;
disp('Total volume in [mm3] of the sample is :');
disp(total_sample_v)
disp('The particle density of the sample [#/mm3] :');
disp(density_sample)

reference_site= 15;
slice_num     = 41;
slice_list    = linspace(0,20,slice_num);
step_r        = slice_list(2);
pair_function = zeros(slice_num,1);
distance_r    = [];
for i = 1:size_arrays
    distance_r(end+1) = sqrt( (x_array(i)-x_array(reference_site))^2+(y_array(i)-y_array(reference_site))^2+(z_array(i)-z_array(reference_site))^2);
end

count1 = 0;
counted = [];
 
for j = 1:slice_num
    count = 0;
    for i = 1:size_arrays
        if distance_r(i) > slice_list(j) && distance_r(i)< slice_list(j)+step_r
            count = count+1;
            count1 = count1+1;
        end 
    end
    counted(end+1)   = count1/size_arrays;
    pair_function(j) = count/(size_arrays*density_sample);
    pair_function(j) = pair_function(j)/(4*pi*(slice_list(j)^2)*step_r);
    
end
%disp(pair_function);
 
fig_paircorrelation = figure('visible','off');
hold on
plot(slice_list,pair_function,'LineWidth',2)
%plot(slice_list,counted)
hold off
ylabel('Pair correlation function $g(r)$','Interpreter','latex')
xlabel('Distance $r$ from an other site [mm]','Interpreter','latex')
exportgraphics(fig_paircorrelation,strcat(true_name,'/pair_correlation.png'));

%% PCA to get asphericity and orientation vector

vector_asphericity = zeros(size_arrays,3); % will be (x,y,z)
vector_asphericity_test = zeros(size_arrays,3);
for i = 1:size_arrays
    vector_asphericity(i,1) = PCA1_x(i) ;
    vector_asphericity(i,2) = PCA1_y(i) ;
    vector_asphericity(i,3) = PCA1_z(i) ;

    vector_asphericity_test(i,1) = x_array(i)+PCA1_x(i) ;
    vector_asphericity_test(i,2) = y_array(i)+PCA1_y(i) ;
    vector_asphericity_test(i,3) = z_array(i)+PCA1_z(i) ;
end

direction_asphericity = figure;
hold on
for i = 1:size_arrays
    %plot3( [x_array(i), vector_asphericity_test(i,1)], [y_array(i), vector_asphericity_test(i,2)], [z_array(i), vector_asphericity_test(i,3)] )
    quiver3(x_array(i), y_array(i), z_array(i), vector_asphericity(i,1),vector_asphericity(i,2),vector_asphericity(i,3),'Linewidth',2 ,'MaxHeadSize',2)
end
hold off
xlabel('$x$ position [mm]','Interpreter','Latex');
ylabel('$y$ position [mm]','Interpreter','Latex');
zlabel('$z$ position [mm]','Interpreter','Latex');

%% Histogram of z projection of asphericity

histogram_asphericity_z_project = figure;
hist(PCA1_z,21)





