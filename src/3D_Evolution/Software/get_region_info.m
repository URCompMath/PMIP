function Omega_new = get_region_info(Gamma,Omega,config)

Omega_new = Omega; 

Nx = size(Omega_new.A_info,1); 
Ny = size(Omega_new.A_info,2); 
Nz = size(Omega_new.A_info,3); 

% Collect pixels and index of a small tube around surfaces
[pixel_set,n_pixels,index] = collect_pixels(Gamma,config);
fprintf('n_pixels = %d\n', n_pixels); 

cc = [-config.image.sizes(1), -config.image.sizes(2), -config.image.sizes(3)];
if(config.image.flag >= 6) % Medical images
     cc = [0 0 0];
end

% Adapt Omega.A_info locally around Gamma_k
if(config.image.flag <= 6)
    for j=1:n_pixels
        lx = pixel_set(j,1);
        ly = pixel_set(j,2);
        lz = pixel_set(j,3);
        k  = pixel_set(j,4);
        
        
        if(1<=lx && lx<=Nx && 1<=ly && ly<=Ny && 1<=lz && lz<=Nz)
            switch index(lx,ly,lz)
                case 1 % = direction + nu_k
                    Omega_new.A_info(lx,ly,lz) = Gamma.orient(k,2);
                case 2 % = direction - nu_k
                    Omega_new.A_info(lx,ly,lz) = Gamma.orient(k,1);
            end
            
            % If index has changed from k_old to k_new adapt I_info and
            % n_info
            if(Omega.A_info(lx,ly,lz) ~= Omega_new.A_info(lx,ly,lz))
                k_old = Omega.A_info(lx,ly,lz);
                k_new = Omega_new.A_info(lx,ly,lz);
                
                
                p = -cc + config.image.h * [lx,ly,lz];
                I0 = func_I(p',config.image.flag,Omega);
                
                Omega_new.I_info(k_new,:) = Omega_new.I_info(k_new,:) + I0';
                Omega_new.I_info(k_old,:) = Omega_new.I_info(k_old,:) - I0';
                Omega_new.n_info(k_new,1) = Omega_new.n_info(k_new,1) + 1;
                Omega_new.n_info(k_old,1) = Omega_new.n_info(k_old,1) - 1;
            end
            
        end
        
        
    end
else
    % Special routing for CT iamge using Omega.Omega_info
    for j=1:n_pixels
        lx = pixel_set(j,1);
        ly = pixel_set(j,2);
        lz = pixel_set(j,3);
        k  = pixel_set(j,4);
        
        if(Omega.Omega_info(Ny+1-ly,lx,lz)==0)
            % Set A_info to -1 outside Omega!
            Omega_new.A_info(lx,ly,lz) = -1;
        else
            % Consider direction only inside Omega!
            if(1<=lx && lx<=Nx && 1<=ly && ly<=Ny && 1<=lz && lz<=Nz)
                switch index(lx,ly,lz)
                    case 1 % = direction + nu_k
                        Omega_new.A_info(lx,ly,lz) = Gamma.orient(k,2);
                    case 2 % = direction - nu_k
                        Omega_new.A_info(lx,ly,lz) = Gamma.orient(k,1);
                end
                
                % If index has changed from k_old to k_new adapt I_info and
                % n_info
                if(Omega.A_info(lx,ly,lz) ~= Omega_new.A_info(lx,ly,lz))
                    k_old = Omega.A_info(lx,ly,lz);
                    k_new = Omega_new.A_info(lx,ly,lz);
                    
                    
                    p = -cc + config.image.h * [lx,ly,lz];
                    I0 = func_I(p',config.image.flag,Omega);
                    
                    Omega_new.I_info(k_new,:) = Omega_new.I_info(k_new,:) + I0';
                    Omega_new.I_info(k_old,:) = Omega_new.I_info(k_old,:) - I0';
                    Omega_new.n_info(k_new,1) = Omega_new.n_info(k_new,1) + 1;
                    Omega_new.n_info(k_old,1) = Omega_new.n_info(k_old,1) - 1;
                end
                
            end
        end
        
        
    end
end

      
for i=1:Omega_new.nr_regions
    if(Omega_new.n_info(i,1)>0)
        Imean = Omega_new.I_info(i,:)/Omega_new.n_info(i,1); 
    else
        Imean = zeros(size(Omega_new.I_info(i,:)));
    end
    Omega_new.coeffs(i,:) = Imean;
end


end

