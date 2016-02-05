function Gamma_new = compute_omega(Gamma)
J = size(Gamma.X,1); 
J0 = size(Gamma.simplices,1); 

Gamma.omega = zeros(J,3);
area_Lambda = zeros(J,1); 

for i=1:J0
    nodes = Gamma.simplices{i,1}.nodes;
    X1 = Gamma.X(nodes(1),:);
    X2 = Gamma.X(nodes(2),:);
    X3 = Gamma.X(nodes(3),:);
    
    % compute nu_sigma, area_sigma
    cross_j    = cross(X1-X3,X2-X3); 
    nu_sigma   = cross_j/norm(cross_j);
    area_sigma = 0.5*norm(cross_j);
    
    
    for k=1:3
        % Update omega
        Gamma.omega(nodes(k),:) = Gamma.omega(nodes(k),:) + area_sigma*nu_sigma;
        area_Lambda(nodes(k),1) = area_Lambda(nodes(k),1) + area_sigma; 
    end
    
end

for k=1:J
    Gamma.omega(k,:) = Gamma.omega(k,:) / area_Lambda(k);
    

end
Gamma_new = Gamma;
end