function Gamma_new = compute_normal_and_area_after_mesh_operation(Gamma)

J0 = size(Gamma.simplices,1); 
J  = size(Gamma.X,1); 

area_Lambda = zeros(J,1);
nu_sigma    = zeros(J0,3); 
area_sigma  = zeros(J0,1);
nabla_sigma = zeros(J0,9); 

omega = zeros(J,3); 
tau1 = zeros(J,3); 
tau2 = zeros(J,3); 
v = zeros(J,3); 


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


% Generate and update output
Gamma_new = Gamma; 
Gamma_new.nu_sigma = nu_sigma;
Gamma_new.area_sigma = area_sigma; 
Gamma_new.area_Lambda = area_Lambda; 
Gamma_new.nabla_sigma = nabla_sigma; 
Gamma_new.omega = omega;
Gamma_new.tau1 = tau1;
Gamma_new.tau2 = tau2; 
Gamma_new.v = v;

end