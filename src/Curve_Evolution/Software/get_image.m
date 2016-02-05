function [Image_new,Surface_new] = get_image(config,Image,Surface)

%% Image data
switch config.dimension
    case 2   % 2D Image
        fprintf('Load image data\n'); 
        switch Image.flag
            case 1  % ABC
                load '../Input/2d/ABC.mat';                                  
            case 2 % Multiphase Gray (3 objects)
                load '../Input/2d/3object.mat';                         
            case 3 % Triple Junctions and Color
                load '../Input/2d/triple.mat';                                
            case 4 % Medical Image
                load '../Input/2d/medical.mat';                                
            case 5 % Flowers, Berkley Database
                load '../Input/2d/flowers.mat';                               
            case 6 % Satellite tracking
                k = config.trackingnr;                                    
                k_str = num2str(k); 
                filename = ['../Input/2d_tracking/image_', k_str, '.mat']; 
                load(filename); 
        end
        
             
    case 3  % Image on a triangulated surface
        fprintf('Load image and surface data\n');

        switch Surface.flag
            case 1 % bunny with three circles
                load '../Input/images_on_surfaces/bunny_balls_noise.mat';  
            case 2 % face001 of Basel Face Database
                load '../Input/images_on_surfaces/face001_color.mat'
            case 3 % face052 of Basel Face Database
                load '../Input/images_on_surfaces/face052_color.mat'
            case 4 % face053 of Basel Face Database
                load '../Input/images_on_surfaces/face053_color.mat'
            case 5 % Earth - longwave radiation
                load '../Input/images_on_surfaces/longwave_radiation.mat'
            case 6 % Earth - net radiation
                load '../Input/images_on_surfaces/net_radiation.mat'
            case 7  % Torus Objects
                load '../Input/images_on_surfaces/torus_obj.mat'

        end
        Surface.tri = tri;
        Surface.nodes = nodes;
        Surface.neigh = neigh;
        I = ImageData;
        
        a_max = ceil(max(nodes))+[1,1,1];
        a_min = floor(min(nodes))-[1,1,1];
        
        Surface.aminmax = [a_min(1) a_max(1) a_min(2) a_max(2) a_min(3) a_max(3)];
        clear tri nodes ImageData
        
        
end


Image.data = I;
Image.draw = I;
Image.approx = zeros(size(I)); 




%% Colormap for gray images
cmap = zeros(256,3);
for i=0:255
    cmap(i+1,:)=[i/255,i/255,i/255];
end
Image.cmap = cmap;


%% Output
Image_new = Image;
Surface_new = Surface;

end