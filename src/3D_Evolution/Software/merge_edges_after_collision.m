function Gamma_new = merge_edges_after_collision(Gamma, assign_simp)

Nnodes = size(Gamma.X,1); 
Nsimp = size(Gamma.simplices,1); 
Nedges = size(Gamma.edges,1); 

mark_edges = zeros(Nedges,1); 
mark_nodes = zeros(Nnodes,1); 

ind_nr = [2,3; 3,1; 1,2]; 

surface_index = zeros(2,1); 
surface_index(1,1) = Gamma.simplices{assign_simp(1,1),1}.index(2);

%% Changes in the nodes and simplices
Nassign = size(assign_simp,1); 
node_pairs = zeros(Nassign,2); 

for i=1:Nassign
    isimp1 = assign_simp(i,1);
    isimp2 = assign_simp(i,2);
    
    sindex = Gamma.simplices{isimp2,1}.index(2);
    if(sindex ~= surface_index(1,1))
        surface_index(2,1) = sindex;
    end
    
        
    kloc1 = assign_simp(i,5) ;
    kloc2 = assign_simp(i,6) ;
    
    % Get node indices and mark nodes for deletion
    k01 = Gamma.simplices{isimp1,1}.nodes(ind_nr(kloc1,1));
    k02 = Gamma.simplices{isimp1,1}.nodes(ind_nr(kloc1,2));
    
    l01 = Gamma.simplices{isimp2,1}.nodes(ind_nr(kloc2,1)) ;
    l02 = Gamma.simplices{isimp2,1}.nodes(ind_nr(kloc2,2)) ;
    
    node_pairs(i,:) = [k01,l02]; 
    
    
    mark_nodes([l01;l02],1) = 1; 
    
    % Get edge indices and mark edges for deletion
    e1 = Gamma.simplices{isimp1,1}.edges(kloc1) ;
    e2 = Gamma.simplices{isimp2,1}.edges(kloc2) ;
    mark_edges(e2,1) = 1; 
    
    
    % Update simplex neighbor info
    Gamma.simplices{isimp1,1}.neigh(kloc1) = isimp2;
    Gamma.simplices{isimp2,1}.neigh(kloc2) = isimp1; 
    
    % Update edge info simp2
    Gamma.simplices{isimp2,1}.edges(kloc2) = e1; 
    
    % Update node info simp2
    Gamma.simplices{isimp2,1}.nodes(ind_nr(kloc2,1)) = k02; 
    Gamma.simplices{isimp2,1}.nodes(ind_nr(kloc2,2)) = k01; 
    
end

%% Compute new value for the nodes
% node_pairs
for i=1:Nassign
    i1 = node_pairs(i,1); 
    i2 = node_pairs(i,2); 
    Gamma.X(i1,:) = (Gamma.X(i1,:) + Gamma.X(i2,:))/2; 
end

%% Update node info other simplices nodes == nodes to be deleted
for i=1:Nsimp
    for k=1:3
        found = 0; 
        l = 1; 
        while(l <= Nassign && found==0)            
            if(Gamma.simplices{i,1}.nodes(k) == node_pairs(l,2))
                Gamma.simplices{i,1}.nodes(k) = node_pairs(l,1);
                found = 1;
            end
            l=l+1;
        end
    end
end

for i=1:Nedges
    for k=1:2
        found = 0; 
        l = 1;
        while(l <= Nassign && found==0)    
            if(Gamma.edges(i,k) == node_pairs(l,2))
                Gamma.edges(i,k) = node_pairs(l,1); 
            end
            l=l+1;
        end
    end
end



%% Delete marked nodes and edges and adapt/reduce other indices
[edges_del,help]=find(mark_edges==1);
[nodes_del,help]=find(mark_nodes==1);


for i=1:Nsimp
    nodes = Gamma.simplices{i,1}.nodes;
    edges0 = Gamma.simplices{i,1}.edges;
    
    for j=1:3
        Gamma.simplices{i,1}.nodes(j) = nodes(j) - size(find(nodes_del<nodes(j)),1);
        Gamma.simplices{i,1}.edges(j) = edges0(j) - size(find(edges_del<edges0(j)),1);
    end
end

for i=1:Nedges
    for j=1:2
        Gamma.edges(i,j) = Gamma.edges(i,j) - size(find(nodes_del < Gamma.edges(i,j)),1);
    end
end

Gamma.X(nodes_del,:) = []; 
Gamma.edges(edges_del,:) = []; 


%% Update surface index 
if(surface_index(2,1) > 0)
    s1 = min(surface_index); 
    s2 = max(surface_index); 
    
    for i=1:size(Gamma.simplices,1)
        if(Gamma.simplices{i,1}.index(2) == s2)
            Gamma.simplices{i,1}.index(2) = s1; 
        end
    end
    Gamma.nr_surfaces = Gamma.nr_surfaces - 1; 
end

%% Output
Gamma_new = Gamma; 




end