function [Gamma_new, mark, mark_delete] = compute_normal_and_area(Gamma,config)


J0 = size(Gamma.simplices,1); 
J  = size(Gamma.X,1); 

area_Lambda = zeros(J,1);
nu_sigma    = zeros(J0,3); 
area_sigma  = zeros(J0,1);
nabla_sigma = zeros(J0,9); 
max_h_sigma = zeros(J0,1); 

omega = zeros(J,3); 
tau1 = zeros(J,3); 
tau2 = zeros(J,3); 
v = zeros(J,3); 

mark = zeros(J0,1);
mark_delete = zeros(J0,2); 

Nmark = 0; 
Nmark_area = 0; 
Nmark_angle = 0; 

Nmark_delete = 0; 
Nmark_delete_area = 0; 
Nmark_delete_angle = 0; 

area_avg = 0; 

phi_max = zeros(J0,1); 
phi_min = zeros(J0,1); 

% Main loop over all simplices, compute nu and area
for j=1:J0
    % compute normal and area of element sigma_j
    nodes = Gamma.simplices{j,1}.nodes;
    X1 = Gamma.X(nodes(1),:);
    X2 = Gamma.X(nodes(2),:);
    X3 = Gamma.X(nodes(3),:); 
        
    % compute nu_sigma, area_sigma
    cross_j         = cross(X1-X3,X2-X3); 
    nu_sigma(j,:)   = cross_j/norm(cross_j);
    area_sigma(j)   = 0.5*norm(cross_j);
    max_h_sigma(j)  = max([norm(X1-X2), norm(X2-X3), norm(X3-X1)]); 
    
    % Add contribution of sigma_j to area_Lambda_k and to omega_k
    for k=1:3
        area_Lambda(nodes(k)) = area_Lambda(nodes(k)) + area_sigma(j);
        omega(nodes(k),:) = omega(nodes(k),:) + area_sigma(j) * nu_sigma(j,:); 
    end
    
    % Compute nabla_sigma
    % Calculate \nabla_s \chi_k for k=1,2,3
    u = X2-X1;
    v = X3-X1;
    nabla_sigma(j,1:3)= -u/(norm(u)^2)-v/(norm(v)^2);
    
    u = X1-X2;
    v = X3-X2;
    nabla_sigma(j,4:6)= -u/(norm(u)^2)-v/(norm(v)^2);

    u = X1-X3;
    v = X2-X3;
    nabla_sigma(j,7:9)= -u/(norm(u)^2)-v/(norm(v)^2);
    
    % Update area_avg
    area_avg = area_avg + area_sigma(j); 
    
    % Compute biggest angle
    v1 = X1-X2;
    v2 = X2-X3;
    v3 = X3-X1;
    a1 = acos(-v1/norm(v1) * (v3/norm(v3))');
    a2 = acos(-v2/norm(v2) * (v1/norm(v1))');
    a3 = acos(-v3/norm(v3) * (v2/norm(v2))');
    
    [a,i3] = max([a1,a2,a3]);
    phi_max(j) = a;
    
    [a,i3] = min([a1,a2,a3]);
    phi_min(j) = a;
            
            
    
end

% Compute omega, tau_1, tau_2
for k=1:J
    omega(k,:) = omega(k,:) / area_Lambda(k); 
    
    v(k,:) = omega(k,:) / norm(omega(k,:)); 
    [val,ind] = min(abs(v(k,:))); 
    
    % v1 unit vector in direction "ind"
    v1 = zeros(1,3); 
    v1(ind) = 1; 
    
    % Adapt v1 such that v1 and v are orthogonal
    v1 = v1 - (v1*v(k,:)')*v(k,:);
    v1 = v1/norm(v1); 
    
    % Get v2
    v2 = cross(v(k,:),v1); 
    
    tau1(k,:) = v1;
    tau2(k,:) = v2; 
    
end

area_avg = area_avg / J0; 

% Compute refinement factor
err_rel = (area_avg - config.mesh.area_desired)/config.mesh.area_desired; 
config.mesh.refine_factor = config.mesh.refine_factor * (1 - err_rel * 0.2);

factor = config.mesh.refine_factor;
factor = max([config.mesh.factor_min, factor]); 
fprintf('Use refinement factor %3.5f, area avg = %3.5f, area avg desired = %3.5f\n', factor, area_avg, config.mesh.area_desired); 



% Mark simplices for refinement or deletion
for j=1:J0
    % Mark for refinement if area too large
    if(area_sigma(j) >= factor*config.mesh.area_desired) 
        Nmark = Nmark + 1;
        Nmark_area = Nmark_area + 1; 
        mark(Nmark) = j;
    else
        % Mark for deletion if area is too small
        if(area_sigma(j) < 0.01*config.mesh.area_desired)
            Nmark_delete = Nmark_delete + 1;
            Nmark_delete_area = Nmark_delete_area + 1;
            mark_delete(Nmark_delete,:) = [j, 1];  % 1: for case "area"
        else
            % Mark for refinement if angle too large
            if(phi_max(j) >= config.mesh.phi_max)
                Nmark = Nmark + 1;
                Nmark_angle = Nmark_angle + 1;
                mark(Nmark) = j;
            else
                % Mark for deletion (edge collapse) if angle too small
                if(phi_min(j) < config.mesh.phi_min)
                    Nmark_delete = Nmark_delete + 1;
                    Nmark_delete_angle = Nmark_delete_angle + 1;
                    mark_delete(Nmark_delete,:) = [j, 2];  %2: for case "angle"
                end
            end
        end
        
    end
    
end
mark = mark(1:Nmark,1); 
mark_delete = mark_delete(1:Nmark_delete,:);     
    

% Print information about refinement / deletion
fprintf('%d simplices marked for refinement.\n', Nmark); 
fprintf('%d simplices have an area >= %2.3f x %2.3f = %2.3f.\n', Nmark_area, factor, config.mesh.area_desired, factor*config.mesh.area_desired); 
fprintf('%d simplices have an angle >= %2.3f deg.\n', Nmark_angle, config.mesh.phi_max * 180/pi); 

fprintf('%d simplices marked for deletion\n', Nmark_delete); 
fprintf('%d simplices have an area < %2.3f x %2.3f = %2.3f.\n', Nmark_delete_area, 0.01, config.mesh.area_desired, 0.01*config.mesh.area_desired); 
fprintf('%d simplices have an angle < %2.3f deg.\n\n', Nmark_delete_angle, config.mesh.phi_min * 180/pi); 




% Generate and update output
Gamma_new = Gamma; 
Gamma_new.nu_sigma = nu_sigma;
Gamma_new.area_sigma = area_sigma; 
Gamma_new.area_Lambda = area_Lambda; 
Gamma_new.nabla_sigma = nabla_sigma; 
Gamma_new.max_h_sigma = max_h_sigma;
Gamma_new.omega = omega;
Gamma_new.tau1 = tau1;
Gamma_new.tau2 = tau2; 
Gamma_new.v = v;

end