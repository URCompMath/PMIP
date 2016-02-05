function Omega = get_initial_region_info(Gamma,config,initial_info)

% Number of regions
Omega.nr_regions = Gamma.nr_surfaces + 1;

% Initial information
c = initial_info.data(:,1:3);
r = initial_info.data(:,4:6);


% Initialize output
Omega.I_info = zeros(Omega.nr_regions,1);
Omega.n_info = zeros(Omega.nr_regions,1);
Omega.coeffs = zeros(Omega.nr_regions,1);

%% Set A_info
Nx = 2*ceil(config.image.sizes(1)/config.image.h);
Ny = 2*ceil(config.image.sizes(2)/config.image.h);
Nz = 2*ceil(config.image.sizes(3)/config.image.h);

cc = [-config.image.sizes(1), -config.image.sizes(2), -config.image.sizes(3)];


% Load region info data (CT examples)    
if(config.image.flag == 6) % Lung TCIA
    cc = [0 0 0];
    Nx = config.image.sizes(1);
    Ny = config.image.sizes(2);
    Nz = config.image.sizes(3);
    
    load('lung_tcia.mat');
    Omega.data = I;
else
    if(config.image.flag == 7) % CT UKR Lung+Heart Segm
        cc = [0 0 0];
        Nx = config.image.sizes(1);
        Ny = config.image.sizes(2);
        Nz = config.image.sizes(3);
        
        load('ct_cd0_part1_01.mat');
        Omega.data = I;
        Omega.Omega_info = Omega_info;
    else
        
        if(config.image.flag == 8) % CT UKR Abdominal Region
            cc = [0 0 0];
            Nx = config.image.sizes(1);
            Ny = config.image.sizes(2);
            Nz = config.image.sizes(3);
            
            load('ct_cd1_part4_03.mat');
            Omega.data = I;
            Omega.Omega_info = Omega_info;
            
        else
            if(config.image.flag == 9) % CT UKR Lung Segm. Splitting
                cc = [0 0 0];
                Nx = config.image.sizes(1);
                Ny = config.image.sizes(2);
                Nz = config.image.sizes(3);
                
                load('ct_cd4_part2_01.mat');
                Omega.data = I;
                Omega.Omega_info = Omega_info;
            end
        end
    end
end
Omega.A_info = zeros(Nx,Ny,Nz);


if(config.image.flag <= 6) 
    for k=1:Nx
        fprintf('%d/%d\n', k,Nx);
        for l=1:Ny
            for m=1:Nz
                p = cc+[k,l,m]*config.image.h;  % point in 3D corresponding to (k,l,m)
                I0 = func_I(p',config.image.flag,Omega);
                
                % Get phase index where (k,l,m) belongs to
                found = 0;
                i=1;
                imain = 1;
                while(i <= Gamma.nr_surfaces && found==0)
                    if(config.multiphase) % Multiphase Adaptation
                        imain = i;
                    end
                    
                    q = p - c(i,:);
                    
                    inside = 0;
                    switch initial_info.flag
                        case 1  % Sphere / Ellipsoid
                            if((q(1)/r(i,1))^2 + (q(2)/r(i,2))^2 + (q(3)/r(i,3))^2 <= 1)
                                inside = 1;
                            end
                            
                        case 2  % Quader / Cube
                            if(-r(i,1)<=q(1) && q(1)<=r(i,1) && -r(i,2)<=q(2) && q(2)<=r(i,2) && -r(i,3)<=q(3) && q(3)<=r(i,3))
                                inside = 1;
                            end
                            
                        case 3 % 2 balls (at +- c )
                            if((q(1)/r(i,1))^2 + (q(2)/r(i,2))^2 + (q(3)/r(i,3))^2 <= 1)
                                inside = 1;
                            end
                        case 4 % torus r(i,3) unused, R=r(i,1), r=r(i,2)
                            if((sqrt(q(1)^2+q(2)^2)-r(i,1))^2+q(3)^2<=r(i,2)^2)
                                inside = 1;
                            end
                        case 5 % cylinder, r(i,3) ununsed, r = r(i,1), h=r(i,2)
                            if(q(1)>= -r(i,2)/2 && q(1)<= r(i,2)/2)
                                if(q(2)^2 + q(3)^2 <= r(i,1)^2)
                                    inside = 1;
                                end
                            end
                            
                            
                    end
                    if(inside == 1)
                        % Get phase index of phase in the interior of Gamma_i
                        phase_index = Gamma.orient(imain,1);
                        if(Gamma.orient(imain,1)== Omega.nr_regions)
                            phase_index = Gamma.orient(imain,2);
                        end
                        found = 1;
                    end
                    i=i+1;
                end
                
                % If found == 0 (not in the interior of any surface, set phase
                % index to nr_regions
                if(found == 0)
                    phase_index = Omega.nr_regions;
                end
                
                % Set A_info, I_info, n_info
                Omega.A_info(k,l,m) = phase_index;
                Omega.I_info(phase_index,:) = Omega.I_info(phase_index,:) + I0';
                Omega.n_info(phase_index,1) = Omega.n_info(phase_index,1) + 1;
                
                
            end
        end
    end
else
    % Special routine for UKR CT images using Omega.Omega_info
    for k=1:Nx
        fprintf('%d/%d\n', k,Nx);
        for l=1:Ny
            for m=1:Nz
                p = [k,l,m];  % point in 3D corresponding to (k,l,m)
                
                if(Omega.Omega_info(Ny+1-round(p(2)),round(p(1)),round(p(3))) == 0)
                    Omega.A_info(k,l,m) = -1;
                else
                    I0 = func_I(p',config.image.flag,Omega);
                    
                    % Get phase index where (k,l,m) belongs to
                    found = 0;
                    i=1;
                    imain = 1;
                    while(i <= Gamma.nr_surfaces && found==0)
                        if(config.multiphase) % Multiphase Adaptation
                            imain = i;
                        end
                        
                        q = p - c(i,:);
                        
                        inside = 0;
                        switch initial_info.flag
                            case 1  % Sphere / Ellipsoid
                                if((q(1)/r(i,1))^2 + (q(2)/r(i,2))^2 + (q(3)/r(i,3))^2 <= 1)
                                    inside = 1;
                                end
                                
                            case 2  % Quader / Cube
                                if(-r(i,1)<=q(1) && q(1)<=r(i,1) && -r(i,2)<=q(2) && q(2)<=r(i,2) && -r(i,3)<=q(3) && q(3)<=r(i,3))
                                    inside = 1;
                                end
                                
                            case 3 % 2 balls (at +- c )
                                if((q(1)/r(i,1))^2 + (q(2)/r(i,2))^2 + (q(3)/r(i,3))^2 <= 1)
                                    inside = 1;
                                end
                            case 4 % torus r(i,3) unused, R=r(i,1), r=r(i,2)
                                if((sqrt(q(1)^2+q(2)^2)-r(i,1))^2+q(3)^2<=r(i,2)^2)
                                    inside = 1;
                                end
                            case 5 % cylinder, r(i,3) ununsed, r = r(i,1), h=r(i,2)
                                if(q(1)>= -r(i,2)/2 && q(1)<= r(i,2)/2)
                                    if(q(2)^2 + q(3)^2 <= r(i,1)^2)
                                        inside = 1;
                                    end
                                end
                                
                                
                        end
                        if(inside == 1)
                            % Get phase index of phase in the interior of Gamma_i
                            phase_index = Gamma.orient(imain,1);
                            if(Gamma.orient(imain,1)== Omega.nr_regions)
                                phase_index = Gamma.orient(imain,2);
                            end
                            found = 1;
                        end
                        i=i+1;
                    end
                    
                    % If found == 0 (not in the interior of any surface, set phase
                    % index to nr_regions
                    if(found == 0)
                        phase_index = Omega.nr_regions;
                    end
                    
                    % Set A_info, I_info, n_info
                    Omega.A_info(k,l,m) = phase_index;
                    Omega.I_info(phase_index,:) = Omega.I_info(phase_index,:) + I0';
                    Omega.n_info(phase_index,1) = Omega.n_info(phase_index,1) + 1;
                    
                end
            end
        end
    end
    
    
end
fprintf('\n\n');

%% Calc coefficients
for i=1:Omega.nr_regions
    if(Omega.n_info(i,1) > 0)
        Imean = Omega.I_info(i,:)/Omega.n_info(i,1);
    else
        Imean = zeros(size(Omega.I_info(i,:)));
    end
    Omega.coeffs(i,:) = Imean;
end

end