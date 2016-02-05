function Gamma_new = improve_mesh2_general(Gamma,config, max_nr_simp, max_loc)

%% Get simplices belonging to nodes
N = size(Gamma.X,1);
nodes_simp = zeros(N,10);
nodes_local_index = zeros(N,10);
nodes_simp_nr = zeros(N,1);

nodes_neigh_vertices = zeros(N,20);
nodes_neigh_nr = zeros(N,1);

Nsimp = size(Gamma.simplices,1);
for i=1:Nsimp
    for j=1:3
        j0 = Gamma.simplices{i,1}.nodes(j);
        nodes_simp_nr(j0,1) = nodes_simp_nr(j0,1) + 1;
        nodes_simp(j0,nodes_simp_nr(j0,1)) = i;
        nodes_local_index(j0,nodes_simp_nr(j0,1)) = j;
        
        j1 = Gamma.simplices{i,1}.nodes(mod(j,3)+1);
        j2 = Gamma.simplices{i,1}.nodes(mod(j+1,3)+1);
        
        nodes_neigh_nr(j0,1) = nodes_neigh_nr(j0,1) + 2;
        nodes_neigh_vertices(j0,nodes_neigh_nr(j0,1)-1) = j1;
        nodes_neigh_vertices(j0,nodes_neigh_nr(j0,1))   = j2;
        
    end
end

%% Improve mesh quality if more than max_nr_simp simplices
index = find(nodes_simp_nr > max_nr_simp);
data = [index, nodes_simp_nr(index,1)];

[nr,ind] = sort(data(:,2),'descend'); 
index = index(ind,:);

n = min([size(index,1), max_loc]); 
index = index(1:n,:); 

fprintf('Improve mesh at %d nodes\n', size(index,1)); 
for i=1:size(index,1)
    j0 = index(i,1);
    
    if(j0 > 0)
        
        kmax = 20;
        k=1;
        n_stop = 7;
        n_remaining = 2*n_stop;
        
        while(k<= kmax && n_remaining > n_stop)
            [Gamma,n_remaining] = improve_mesh2(Gamma,j0); 
            
            k=k+1;
        end
        index_neigh = nodes_neigh_vertices(j0,1:nodes_neigh_nr(j0,1));
        for k=1:size(index_neigh,2)
            k0 = index_neigh(k);
            k1 = find(index==k0);
            index(k1,1) = 0;
        end
        
    end
    
end
Gamma = compute_normal_and_area(Gamma,config);
Gamma_new = Gamma;

end


