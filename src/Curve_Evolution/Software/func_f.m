function y = func_f(i_node,Gamma,Image,Omega,config,curve_nr)
if(config.dimension==2)
    fprintf('Use funcf (c-code) for dim=2\n');
    pause
end

i_simp = Gamma.closest_simp(i_node,1);
I0 = Image.data(i_simp,:)';


% Region based Active Contours Chan Vese with possible Multiphase and Triple Junctions
ci_1 = Omega.coeffs(Gamma.orient(curve_nr,1),:);
ci_2 = Omega.coeffs(Gamma.orient(curve_nr,2),:);

r1 = Gamma.orient(curve_nr,1);
r2 = Gamma.orient(curve_nr,2);

% lambda: columns --> correspond to color components
% rows--------------> correspond to region
lambda = config.method.lambda;




if( size(Image.data,2)==3)  % colored image
    switch config.method.color
        case 1 % RGB
            % normalized differences
            v1 = (ci_1 - I0') / 255;
            v2 = (ci_2 - I0') / 255;
            
            % Compute Output
            y = 0;
            for k=1:3
                y = y + (lambda(r2,k) * v2(k)*v2(k) - lambda(r1,k) * v1(k)*v1(k) ) ;
            end
            
        case 2 % CB
            % lambda = [lambda_c, lambda_b]
            
            % Transform I0 to cb
            Icb = my_rgb2cb(I0);
            
            % Differences
            v1c = ci_1(1:3) - Icb(1:3)';
            v1b = ci_1(4) - Icb(4);
            v2c = ci_2(1:3) - Icb(1:3)';
            v2b = ci_2(4) - Icb(4);
            
            % Compute Output
            y =  ( lambda(r2,1) * (v2c*v2c') - lambda(r1,1) * (v1c*v1c')) +  (lambda(r2,2) * v2b^2 - lambda(r1,2) * v1b^2);
            
        case 3 % HSV
            % lambda = [lambda_h, lambda_s, lambda_v]
            
            
            % Transform I0 to hsv
            Ihsv = my_rgb2hsv(I0);
            % convert hue from [0,2pi] to S^1
            Ihsv = [cos(Ihsv(1)); sin(Ihsv(1)); Ihsv(2); Ihsv(3)];
            % Differences
            
            % Hue
            v1h = ci_1(1:2) - Ihsv(1:2)';
            v2h = ci_2(1:2) - Ihsv(1:2)';
            
            % Saturation
            v1s = ci_1(3) - Ihsv(3);
            v2s = ci_2(3) - Ihsv(3);
            
            % Value
            v1v = ci_1(4) - Ihsv(4);
            v2v = ci_2(4) - Ihsv(4);
            
            y =  (lambda(r2,1) * (v2h*v2h') - lambda(r1,1) * (v1h*v1h')) +  (lambda(r2,2) * v2s^2 - lambda(r1,2) * v1s^2) +  (lambda(r2,3) * v2v^2 - lambda(r1,3) * v1v^2);
            
            
            
    end
    
    
    
else
    % scalar image
    v1 = (ci_1-I0)/255;
    v2 = (ci_2-I0)/255;
    y = ( lambda(r2,1) *  v2^2 - lambda(r1,1) * v1^2);
end






end

