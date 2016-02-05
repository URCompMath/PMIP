function Gamma_new = delete_surfaces(Gamma,mark,config)
%-------------------------------------------------------------------------%
% Function:           delete_surfaces
%
% Description:        Delete surfaces which are marked for deletion
%
%-------------------------------------------------------------------------%


N = size(Gamma.simplices,1);
N_nodes = size(Gamma.X,1);
N_edges = size(Gamma.edges,1);

% Initialize
X_new = Gamma.X;
simplices_new = Gamma.simplices;
edges_new = Gamma.edges;
area_info_new = Gamma.area_info;

if(sum(sum(mark))>0)
    % Set info_simplices and info_nodes to 1 if simplex and/or node has to be
    % deleted
    info_simplices = zeros(1,N);
    info_nodes     = zeros(1,N_nodes);
    info_edges     = zeros(1,N_edges);
    
    for i=1:N
        if(mark(1,Gamma.simplices{i,1}.index(2))==1)  %simplex is marked for deletion
            info_simplices(1,i)=1;
            for j=1:3
                info_nodes(1,Gamma.simplices{i,1}.nodes(j))=1;
                info_edges(1,Gamma.simplices{i,1}.edges(j))=1;
            end
        end
    end
    
    % Get index vector of nodes, simplices and surfaces which have to be
    % deleted
    [help,nodes_del]  = find(info_nodes==1);
    [help,sim_del]    = find(info_simplices==1);
    [help,edges_del]  = find(info_edges==1);
    [help,surface_del]= find(mark==1);
    
    
    % Adapt output
    for i=1:N
        index = Gamma.simplices{i,1}.index(2);
        if(mark(1,index)==0)
            nodes = Gamma.simplices{i,1}.nodes;
            neigh = Gamma.simplices{i,1}.neigh;
            edges0 = Gamma.simplices{i,1}.edges;
            
            nr = size(find(surface_del<index),2);
            if(config.multiphase)
                nr = size(find(surface_del<index),1);
            end
            
            
            simplices_new{i,1}.index(2) = index - nr;
            for j=1:3
                simplices_new{i,1}.nodes(j) = nodes(j) - size(find(nodes_del<nodes(j)),2);
                simplices_new{i,1}.neigh(j) = neigh(j) - size(find(sim_del<neigh(j)),2);
                simplices_new{i,1}.edges(j) = edges0(j) - size(find(edges_del<edges0(j)),2);
            end
            
        end
    end
    
    for i=1:N_edges
        if(isempty(find(edges_del==i))==1)
            for j=1:2
                edges_new(i,j)= Gamma.edges(i,j) - size(find(nodes_del<Gamma.edges(i,j)),2);
            end
        end
    end
    
    
    
    % Delete rows in simplices_new and X_new
    simplices_new(sim_del,:)=[];
    X_new(nodes_del,:)=[];
    area_info_new(:,surface_del)=[];
    edges_new(edges_del,:)=[];
end

% Set nr_surfaces_new
nr_surfaces_new = Gamma.nr_surfaces - sum(mark);

% Generate output
Gamma_new = Gamma;
Gamma_new.simplices = simplices_new;
Gamma_new.X = X_new;
Gamma_new.area_info = area_info_new;
Gamma_new.edges = edges_new;
Gamma_new.nr_surfaces = nr_surfaces_new;

end

