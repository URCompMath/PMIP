function Gamma_out = area_compute_and_zero_check(Gamma,X_new,simplices_new,edges_new)
% Compute new simplices area and normals, check for area=0 simplices

% Preparation
Gamma_new = Gamma;
Gamma_new.X = X_new;
Gamma_new.simplices = simplices_new;
Gamma_new.edges = edges_new;

% Compute new area and normals
Gamma_new = compute_normal_and_area_after_mesh_operation(Gamma_new); % version without checking for refinement

% Area == 0?
i=1;
found = 0;
while(found == 0 && i<= size(Gamma_new.area_sigma,1))
    if(abs(Gamma_new.area_sigma(i,1)) < 1e-10)
        found = 1;
    end
    i=i+1;
end

if(found==0)
    Gamma_out = Gamma_new;
else
    Gamma_out = Gamma;
    fprintf('Found a simplex with area = 0, reset refine/coarsening, use old surface\n'); 
end

end