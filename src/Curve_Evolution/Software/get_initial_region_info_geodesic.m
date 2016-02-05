function Omega  = get_initial_region_info_geodesic(Gamma,initial_information,Surface,Image,config)
%% Preparations
dim = size(Gamma.X,2);

% Initial information
c = initial_information.center;
r = initial_information.radius;

% nr of main curves
Nmain = size(r,3);

% color image flag
if(dim==2)
    if(size(size(Image.data),2)==3)  % Nx x Ny (x Nz)
        color_flag = 1;
    else
        color_flag = 0;
    end
else
    if(size(Image.data,2)==3)   % Nsimp x 1, Nsimp x 3
        color_flag = 1;
    else
        color_flag = 0;
    end
end

% Initialize output
NR = Nmain+1;
if(color_flag)
    switch config.method.color
        case 1 % RGB, 3 elements: red,green,blue
            Omega.I_info = zeros(NR,3);
            Omega.coeffs = zeros(NR,3); 
        case 2 % CB, 4 elements: chrom1, chrom2, chrom3, brightness
            Omega.I_info = zeros(NR,4);
            Omega.coeffs = zeros(NR,4); 
        case 3 % HSV, 4  elements: 2 for hue (consider S^2 instead of [0,2pi], 
               % 1 for saturation, 1 for value
            Omega.I_info = zeros(NR,4);
            Omega.coeffs = zeros(NR,4); 
        case 4 % Lab, 3 elements: L,a,b
            Omega.I_info = zeros(NR,3);
            Omega.coeffs = zeros(NR,3);            
    end
    
else
    Omega.I_info = zeros(NR,1); % scalar (gray-value) image
    Omega.coeffs = zeros(NR,1);
end
Omega.n_info = zeros(NR,1);

Omega.A_info = zeros(size(Image.data,1),1); % Image is stored as N_simp x 1 matrix!!! 

%% Set A_info, I_info, n_info for dim == 3 (Cases are distinguished already in main.m) 

% Set A_info, I_info, n_info
Omega = get_region_info_geodesic(Gamma, Surface, Omega, Image, config, size(Surface.tri,1), 1);




%% Calc coefficients

for i=1:NR
    % Imean is the mean of the image (scalar and rgb) or a "mean" of hsv or
    % cb values, 
    if(Omega.n_info(i,1)>0)
        Imean = Omega.I_info(i,:)/Omega.n_info(i,1);
    else
        Imean = zeros(size(Omega.I_info(i,:))); 
    end
    
    if(color_flag)
        switch config.method.color 
            case 1 % RGB
                Omega.coeffs(i,:) = Imean; 
            case 2 % CB
                
                % Imean(4) stores brightness mean 
                Omega.coeffs(i,4) = Imean(4); 
                
                % Imean(1:3) is the right direction but does not lie on S^2
                % needs to be normalized
                if(norm(Imean(1:3))>0)
                    Omega.coeffs(i,1:3) = Imean(1:3)/norm(Imean(1:3));
                end
                
            case 3 % HSV
                if(norm(Imean(1:2))>0)
                    Omega.coeffs(i,1:2) = Imean(1:2)/norm(Imean(1:2)); % hue: project I(1:2) to S^1
                else
                    Omega.coeffs(i,1:2) = [0,0]; 
                end
                Omega.coeffs(i,3) = Imean(3);   % saturation
                Omega.coeffs(i,4) = Imean(4);   % value
                                
                
        end
                
    else
        Omega.coeffs(i,:) = Imean;
    end
end

end

