function [Eint,Eext] = compute_energies(config,Gamma,Image,Omega,Surface)
% 2D only, cannot be used for images on surfaces (see line 27ff)

% Eint = sigma * sum of lengths
Eint = 0;

kmain = size(Gamma.index_info,1);
for k=1:kmain
    for l=1:Gamma.nr_curves(k,1)
        Eint = Eint + Gamma.length(k,l);
    end
end
Eint = config.method.sigma * Eint;

% Eext = int_Omega (u-u0)^2 dx
Eext = 0;

step = config.method.sigma_step; 

lambda = config.method.lambda;
if(size(size(Image.data),2)==3)  % color image
    switch config.method.color
        case 1 % RGB
            % compute contribution of each pixel to Eext
            for k=1:step:Image.sizes(2)
                for l=1:step:Image.sizes(1)                    
                    z = [l, Image.sizes(2)+1-k];
                    
                    I0 = func_I(z,Image,Surface);
                    
                    % normalized differences
                    v = (Omega.coeffs(Omega.A_info(k,l),:) - I0') / 255;
                    
                    % Compute contribution to Eext
                    for d=1:3
                        Eext = Eext + lambda(Omega.A_info(k,l),d) *  v(d)*v(d) * step^2;
                    end
                end
            end


        case 2 % CB
            % compute contribution of each pixel to Eext
            for k=1:step:Image.sizes(2)
                for l=1:step:Image.sizes(1)
                   z = [l, Image.sizes(2)+1-k];
                   
                   I0 = func_I(z,Image,Surface);

                   % Transform I0 to cb
                   Icb = my_rgb2cb(I0);
                   
                   % Squared Differences
                   ci = Omega.coeffs(Omega.A_info(k,l),:);
                   vc = ci(1:3) - Icb(1:3)';
                   vb = ci(4) - Icb(4);
                   
                   % Compute contribution to Eext
                   Eext = Eext + (lambda(Omega.A_info(k,l),1) * (vc*vc') + lambda(Omega.A_info(k,l),2) * vb^2) * step^2;
                    
                end
            end
            
        case 3 % HSV
            % compute contribution of each pixel to Eext
            for k=1:step:Image.sizes(2)
                for l=1:step:Image.sizes(1)
                    z = [l, Image.sizes(2)+1-k];

                    I0 = func_I(z,Image,Surface);

                    % Transform I0 to hsv
                    Ihsv = my_rgb2hsv(I0);
                    
                    % convert hue from [0,2pi] to S^1
                    Ihsv = [cos(Ihsv(1)); sin(Ihsv(1)); Ihsv(2); Ihsv(3)];
                    
                    % Squared Differences
                    ci = Omega.coeffs(Omega.A_info(k,l),:);
                    vh = ci(1:2) - Ihsv(1:2)';   % hue
                    vs = ci(3) - Ihsv(3);        % saturation
                    vv = ci(4) - Ihsv(4);        % value
                    
                    % Compute contribution to Eext
                    Eext = Eext + (lambda(Omega.A_info(k,l),1) * (vh*vh') + lambda(Omega.A_info(k,l),2) * vs^2 + lambda(Omega.A_info(k,l),3) * vv^2) * step^2 ;
                    
                end
            end
            

    end
else
    % scalar image
    % compute contribution of each pixel to Eext
    for k=1:step:Image.sizes(2)
        for l=1:step:Image.sizes(1)
            z = [l, Image.sizes(2)+1-k];

            I0 = func_I(z,Image,Surface);

            % scalar image
            v = (Omega.coeffs(Omega.A_info(k,l),:)-I0)/255;
            
            % Compute contribution to Eext
            Eext = Eext + (lambda(Omega.A_info(k,l),1)* v^2) * step;
        end
    end
            

end
            

end
        
        
        
        
        
        
        
        