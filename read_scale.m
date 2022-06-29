function [out] = read_scale(name_image)
    if nargin==0      
    error('There are not enought input arguments');    
    elseif nargin>1
    error('Too many input arguments');             
    elseif nargin==1 
        picture                         = imread(name_image);
        figure;
        imshow(picture);
        [points_scale_x,points_scale_y] = ginputWhite(2);
        x1                              = floor(points_scale_x(1)); 
        x2                              = floor(points_scale_x(2));
        y1                              = floor(points_scale_y(1));
        y2                              = floor(points_scale_y(2));
        physical_size                   = input('What is the physical size of the scale (in mm)?');
        gs                              = rgb2gray(picture); %Shading of grey
        BW                              = imbinarize(gs,0.25); %Binarization
        BW1                             = BW(y1:y2,x2:x1); %Resizing the binary image with the new coordinates
        dimensions_cache                = [];
        [numRows,numCols]               = size(BW1);
        for i = 1:numRows
            length_cache = 0;
            for j = 1:numCols
                if BW1(i,j) ~=0
                    length_cache = length_cache+1;
                end
            end
            dimensions_cache(end+1)=length_cache;
        end
        disp(strcat('The image is a resolution of  ', num2str(physical_size / max(dimensions_cache)),' mm per px'));
    end
end




