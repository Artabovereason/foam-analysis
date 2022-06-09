close all
clear
clc

%% Import parameters

A = fullfile('/Users/Name/Documents/matlab_param.xlsx');
param = readtable(A);
alpha_800 = param{5,5};
gamma_26 = 0.026;

ext        = '*.tif';
count      = 0;
Nom_essai  = {zeros};
dF_cell={zeros};
l_cell={zeros};
abs_cell = {zeros};
ord_cell = {zeros};
moy_abs = {zeros};
moy_ord = {zeros};
err_abs = {zeros};
err_ord = {zeros};


%% Treatment
while 1
    count=count+1;
    %% CHOOSE FILES
    cd /Users/Name/Documents/;
    [f, rep] = uigetfile('*.tif','Image tomo à traiter :');
    chemin = fullfile(rep, ext);
    list = dir(chemin);
    
    %% RESIZE TOMOGRAPH IMAGES
    image = fullfile(rep, list(1).name);
    I = imread(image); %Original image
    gs = rgb2gray(I); %Shading of grey
    I_BW = imbinarize(gs,0.25); %Binarization
    figure;
    imshow(I_BW);
    
    [xi,yi] = ginputWhite(2); %Using the command 'ginputWhite' to define a ROI
    x1 = floor(xi(1)); % x1, x2, y1, y2 new ROI image coordinates
    x2 = floor(xi(2));
    y1 = floor(yi(1));
    y2 = floor(yi(2));
    
    for n = 1:length(list)
        image = fullfile(rep, list(n).name); % Open the n images one after the other
        I = imread(image); %Original image
        gs = rgb2gray(I); %Shading of grey
        BW = imbinarize(gs,0.25); %Binarization
        
        BW1 = BW(y1:y2,x1:x2); %Resizing the binary image with the new coordinates
        figure;
        imshow(BW1);
       
           
    %% Consider the ribbon as a pixel line corresponding to the centre of the ribbon
    nblign=size(BW1,1);
    xres=zeros(nblign,1);yres=xres; %creation matrix of zeros 
    for j=1:nblign
        [dum,xx]=find(BW1(j,:)==1);
        xres(j)=floor(mean(xx));
        yres(j)=j;
    end
    
    hold on
    plot(xres,yres)
    
        %MEASUREMENTS OF EXTREMUMS [m]
        pixel_m = 0.0000315; % resolution [m]
        max_x_m = max(xres)*pixel_m;
        min_x_m = min(xres)*pixel_m;
        max_y_m = max(yres)*pixel_m;
        min_y_m = min(yres)*pixel_m;
        
        %FIND THE y CORRESPONDING TO max_x AND min_y
        ymin_x=yres(round(mean(find(xres==min(xres)))));
        ymax_x=yres(round(mean(find(xres==max(xres)))));
        ymin_x_m=yres(round(mean(find(xres==min(xres)))))*pixel_m;
        ymax_x_m=yres(round(mean(find(xres==max(xres)))))*pixel_m;
        
        
        %CALCULATION OF dF [m]
        dF_cell{n} = max_x_m - min_x_m;
      
        %CALCULATION OF  ell [m]
        diff= abs((ymax_x_m - ymin_x_m));
        l_cell{n} = diff;
 
            abs_cell{n} = alpha_800./(gamma_26*l_cell{n}.^2); %Calculation abscissa: eta = alpha /(gamma*l^2)
            absc = cell2mat(abs_cell);
            ord_cell{n} = sqrt(3)*dF_cell{n}./l_cell{n}; %Calculation ordinate = sqrt(3)*dF/l
            ord = cell2mat(ord_cell);

        
    end
    op=upper(input('Voulez-vous traiter un autre fichier ? (Y/N)','s'));
    
    [filepath,~,~] = fileparts(rep);
    [~,Nom_essai{count}] = fileparts(filepath);
    

    % CALCULATION OF AVERAGE
    moy_abs{count} = mean(absc);
    moy_ord{count} = mean(ord);
    
    err_abs{count} = std(absc); %error abscissa per 100 cuts 
    err_ord{count} = std(ord); %error ordinate per 100 cuts 
    
    if upper(op) == 'N', break, end
    close all
end

%% SAVE TABLE PARAMETERS

list_nomsVariables = {'Nom_essai','Moyenne_abscisses','Moyenne_ordonnees','erreur_abscisses_100_coupes','erreur_ordonnees_100_coupes'};
T = table(Nom_essai(:),moy_abs(:),moy_ord(:),err_abs(:),err_ord(:),'VariableNames',list_nomsVariables);

writetable(T,'/Users/Name/Documents/.txt');

