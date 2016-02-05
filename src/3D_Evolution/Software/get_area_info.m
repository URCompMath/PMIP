function [area_info,mark] = get_area_info(Gamma,config)

nr_surfaces = min([Gamma.nr_surfaces, config.max_nr_surfaces]);

% Initialize area_info and mark
area_info = zeros(1,config.max_nr_surfaces);
mark      = zeros(1,config.max_nr_surfaces); 

if(config.multiphase)
    area_info = zeros(Gamma.nr_surfaces(1,1),1); 
    mark = zeros(Gamma.nr_surfaces(1,1),1); 
end

N = size(Gamma.simplices,1); 

% Add area of single simplexes to area_info(1,l) for l=simplices{i,1}.index
for i=1:N
    nodes = Gamma.simplices{i,1}.nodes;
    X1 = Gamma.X(nodes(1),:);
    X2 = Gamma.X(nodes(2),:);
    X3 = Gamma.X(nodes(3),:); 
    area_info(Gamma.simplices{i,1}.index(1),Gamma.simplices{i,1}.index(2)) = area_info(Gamma.simplices{i,1}.index(1),Gamma.simplices{i,1}.index(2)) + 0.5*norm(cross(X1-X3,X2-X3));
    
end

for i=1:nr_surfaces
    if(config.multiphase == 0)
        ainfo = area_info(1,i); 
    else
        ainfo = sum(area_info(i,:)); 
    end
    
    if(ainfo<config.delete.area_tol)
        mark(1,i) = 1;
        fprintf('Will delete surface nr. (1,%d). Area = %2.2f < tol = %2.2f\n', i, area_info(1,i), config.delete.area_tol);
    end
end

    

end