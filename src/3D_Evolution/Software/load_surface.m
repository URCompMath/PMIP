function Gamma_new = load_surface(Gamma,config)
% Load surface
load(config.initial_surface_filename);

Gamma.simplices = simplices;
Gamma.X = X; 
Gamma.edges = edges;

% Init orientation
Gamma.orient = zeros(Gamma.nr_surfaces,2); 
for i=1:Gamma.nr_surfaces
    Gamma.orient(i,:) = [i, Gamma.nr_surfaces+1]; 
end

Gamma_new = Gamma; 
end