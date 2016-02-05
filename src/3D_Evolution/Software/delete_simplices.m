function Gamma_new = delete_simplices(Gamma,mark)
%-------------------------------------------------------------------------%
% Function:           delete_simplices
%
% Description:        Delete surfaces which are marked for deletion 
%
%-------------------------------------------------------------------------%
simplices = Gamma.simplices;
X = Gamma.X;
edges = Gamma.edges; 

N = size(simplices,1); 
N_nodes = size(X,1);
N_edges = size(edges,1);

% Set info_simplices,info_nodes, info_edges to 1 if it has to be deleted
info_simplices = mark;
info_nodes     = -1*ones(N_nodes,1);
info_edges     = -1*ones(N_edges,1); 

for i=1:N
    if(mark(i)==1)  % simplex is marked for deletion
        for j=1:3
            if(info_nodes(simplices{i,1}.nodes(j))==-1)
                info_nodes(simplices{i,1}.nodes(j))=1;
            end
            if(info_edges(simplices{i,1}.edges(j))==-1)
                info_edges(simplices{i,1}.edges(j))=1;
            end
        end
    else
        for j=1:3
            info_nodes(simplices{i,1}.nodes(j))=0;
            info_edges(simplices{i,1}.edges(j))=0;
        end
    end
end

% Get index vector of nodes, simplices and surfaces which have to be
% deleted
[nodes_del,help]  = find(info_nodes==1);
[sim_del,help]    = find(info_simplices==1);
[edges_del,help]  = find(info_edges==1);



% Initialize output
X_new = X;
simplices_new = simplices;
edges_new =edges;


% Compute/Adapt output 
for i=1:N
    if(mark(i,1)==0)
        nodes = simplices{i,1}.nodes;
        neigh = simplices{i,1}.neigh;
        edges0 = simplices{i,1}.edges;
        
        % Set neigh to -1 if it is a free edge
        for j=1:3
            
            if(simplices{i,1}.neigh(j)>0 && mark(simplices{i,1}.neigh(j))==1)
                simplices_new{i,1}.neigh(j)=-1;
            else
                if(simplices{i,1}.neigh(j)>0)
                    simplices_new{i,1}.neigh(j) = neigh(j) - size(find(sim_del<neigh(j)),1);
                end
            end
        end  
        
        
        
        for j=1:3
            simplices_new{i,1}.nodes(j) = nodes(j) - size(find(nodes_del<nodes(j)),1);         
            simplices_new{i,1}.edges(j) = edges0(j) - size(find(edges_del<edges0(j)),1);
        end
     
    end
end

        
for i=1:N_edges
    if(isempty(find(edges_del==i))==1)
        for j=1:2
            edges_new(i,j)= edges(i,j) - size(find(nodes_del<edges(i,j)),1);
        end
    end
end
        
% Delete rows in simplices_new and X_new
simplices_new(sim_del,:)=[];
X_new(nodes_del,:)=[];
edges_new(edges_del,:)=[];

% Generate output
Gamma_new = Gamma;
Gamma_new.simplices = simplices_new;
Gamma_new.X = X_new;
Gamma_new.edges = edges_new;


end

        