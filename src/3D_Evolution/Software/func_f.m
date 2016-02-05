function y = func_f(z,Gamma,Omega,config,surface_nr)


switch config.method.pde_flag
    case 0  % Mean Curvature flow
        y = 0;
    case 1  % Image Segmentation (Chan Vese)
        r1 = Gamma.orient(surface_nr,1);
        r2 = Gamma.orient(surface_nr,2);
        
        ci_1 = Omega.coeffs(r1,:) ;
        ci_2 = Omega.coeffs(r2,:) ;
        
        % Currently only gray-scaled images, range = [0,255]
        I0 = func_I(z,config.image.flag,Omega); 

        v1 = (ci_1 - I0) / 255; 
        v2 = (ci_2 - I0) / 255; 
        
        y = -config.method.lambda*( v2^2 - v1^2); 
end


if(config.image.flag == 8)
    sat = 1.0; 
    y = sign(y)*min([sat, abs(y)]); 
end
    
    
end