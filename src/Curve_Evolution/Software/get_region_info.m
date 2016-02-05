function Omega_new = get_region_info(Gamma,Omega,config,Image,Surface)

dim = config.dimension;
Omega_new = Omega;

% colored image flag
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

% Initialize Omega_new.coeffs
NR = size(Omega_new.I_info,1);
if(color_flag)
    if(dim==2)
        I0 = [Image.data(:,:,1), Image.data(:,:,2), Image.data(:,:,3)];
    end
       
    switch config.method.color
        case 1 % RGB
            Omega_new.coeffs = zeros(NR,3); 
        case 2 % CB
            Omega_new.coeffs = zeros(NR,4); 
        case 3 % HSV
            Omega_new.coeffs = zeros(NR,4); 
    end
else
    if(dim==2)
        I0 = Image.data;
    end
    Omega_new.coeffs = zeros(NR,1); % scalar case
    config.method.color = 0;
end


%% Update Omega.A_info, Omega.I_info, Omega.n_info
switch dim
    case 2
        Nmain = size(Gamma.index_info,1);
        for k=1:Nmain
            % Call c-code
            [pixel_set,n_pixels,index] = collectpixels(Gamma.X,Gamma.neigh,Gamma.nr_curves,Gamma.index_info,k,Gamma.mark,config.width_band,Image.sizes);          
            pixel_set = pixel_set(1:n_pixels,:); 
            
            test=regions_2d(pixel_set,index,Omega.A_info,Omega.I_info,Omega.n_info,Gamma.orient,double(I0),config.method.color,k);
            
        end
        
    case 3
        % Set A_info, I_info, n_info
        Omega = get_region_info_geodesic(Gamma, Surface, Omega, Image, config, config.width_band, 0);

        
end

%% Calc coeffs
Omega_new = Omega;

for i=1:NR
    % Imean is the mean of the image (scalar and rgb) or a "mean" of hsv or
    % cb values, 

    if(Omega_new.n_info(i,1) > 0)
        Imean = Omega_new.I_info(i,:)/Omega_new.n_info(i,1);
    else
        Imean = zeros(1,size(Omega_new.I_info,2));
    end
    
    if(color_flag)
        % Colored image
        switch config.method.color
            
            case 1 % RGB
                Omega_new.coeffs(i,:) = Imean;
            case 2 % CB
                % Imean(4) stores the brightness mean
                Omega_new.coeffs(i,4) = Imean(4);
                
                % Imean(1:3) is the right direction but does not lie on S^2
                % needs to be normalized                
                if(norm(Imean(1:3))>0)
                    Omega_new.coeffs(i,1:3) = Imean(1:3)/norm(Imean(1:3));
                end
            case 3 % HSV
                 if(norm(Imean(1:2))>0)
                    Omega_new.coeffs(i,1:2) = Imean(1:2)/norm(Imean(1:2)); % hue: project I(1:2) to S^2
                else
                    Omega_new.coeffs(i,1:2) = [0,0]; 
                end
               
                Omega_new.coeffs(i,3) = Imean(3); % saturation
                Omega_new.coeffs(i,4) = Imean(4); % value
                
        end
        
    else
        % Scalar case
        Omega_new.coeffs(i,:) = Imean;
    end

end
end

