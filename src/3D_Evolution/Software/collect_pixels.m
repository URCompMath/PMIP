function [pixel_set,n_pixels,index] = collect_pixels(Gamma,config)


w = config.width_band;
sizes = config.image.sizes;
cc = -sizes;

J0 = size(Gamma.simplices,1);

Nx = 2*ceil(config.image.sizes(1)/config.image.h);
Ny = 2*ceil(config.image.sizes(2)/config.image.h);
Nz = 2*ceil(config.image.sizes(3)/config.image.h);

if(config.image.flag >= 6) % Medical images
    cc = [0 0 0];
    Nx = config.image.sizes(1);
    Ny = config.image.sizes(2);
    Nz = config.image.sizes(3);
end


% initialize output
n_pixels = 0;
pixel_set = zeros(0,4); %(k,l,m,surface_nr)
index = zeros(Nx,Ny,Nz);


% Collect pixels in square with width 2*w around Gamma.X
for i=1:J0
    nodes = Gamma.simplices{i,1}.nodes;
    nu = Gamma.nu_sigma(i,:);
    
    
    for k=1:3
        for l=1:w
            % in direction -nu  

            y = Gamma.X(nodes(k),:) + l*config.image.h*nu;
            y = y-cc;
            yk = ceil(y(1)/config.image.h);
            yl = ceil(y(2)/config.image.h);
            ym = ceil(y(3)/config.image.h);
            
            inside = 1;
            if(yk<=0 || yk>Nx || yl<=0 || yl>Ny || ym<=0 || ym>Nz)
                inside = 0;
            end

            
            if(inside)
                if(index(yk,yl,ym)==0)
                    index(yk,yl,ym)=1;
                    n_pixels = n_pixels + 1;
                    pixel_set(n_pixels,:) = [yk,yl,ym,Gamma.simplices{i,1}.index(1)];
                end
            end
            % in direction -nu
            y = Gamma.X(nodes(k),:) - l*config.image.h*nu;
            y = y-cc;
            yk = ceil(y(1)/config.image.h);
            yl = ceil(y(2)/config.image.h);
            ym = ceil(y(3)/config.image.h);
            
            inside = 1;
            if(yk<=0 || yk>Nx || yl<=0 || yl>Ny || ym<=0 || ym>Nz)
                inside = 0;
            end
            
            if(inside)
                if(index(yk,yl,ym)==0)
                    index(yk,yl,ym)=2;
                    n_pixels = n_pixels + 1;
                    pixel_set(n_pixels,:) = [yk,yl,ym,Gamma.simplices{i,1}.index(1)];
                end
            end
            
        end
    end
end
end
