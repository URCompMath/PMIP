function [Gamma_new, index_new, index_edges_new, index_edges_corners, index_simp_new] = insert_new_points(Gamma,insert,index,index_edges,index_simp)
% Input
% Gamma                struct, current surface structure
% insert               Ninsert x 5, cols: simp index, edge index to be
%                      divided, node indices (col 3&4), number of nodes 
%                      to be inserted 
% index                N x 1, indices of node indices in the circle of
%                      free edges
% index_edges          N x 1, indices of edges in the circle of free edges
% index_simp           N x 1, indices of simplices in the circle of free
%                      edges

% Output
% Gamma_new            struct, updated surface structure
% index_new            Nnew x 1, indices of node indices in the circle of
%                      free edges, update of index
% index_edges_new      Nnew x 1, number of edge indices in the circle of 
%                      free edges, update of index_edges
% index_edges_corners  Nnew x 2, indices of corner points corresponding to
%                      edge of index_edges
% index_simp_new       Nnew x 1, indices of simplices in the circel of free
%                      edges, update of index_simp

Nsimp  = size(Gamma.simplices,1);
Nnodes = size(Gamma.X,1);
Nedges = size(Gamma.edges,1) ;


Ninsert = size(insert,1); 
ind_nr = [2,3;3,1;1,2]; 




for i=1:Ninsert
    isimp  = insert(i,1) ;
    iedge  = insert(i,2);
    
    % Get local index of free edge iedge
    found = 0; 
    k=1;
    while(k<=3 && found==0)
        if(Gamma.simplices{isimp,1}.edges(k) == iedge)
            found = k;
        end
        k=k+1;
    end
    k0 = found; 
    ik0 = Gamma.simplices{isimp,1}.nodes(k0) ;   
    
    % Get the two other nodes
    inode1 = Gamma.simplices{isimp,1}.nodes(ind_nr(k0,1)) ;
    inode2 = Gamma.simplices{isimp,1}.nodes(ind_nr(k0,2)) ;
    
    
    
    % Update number of simplices, nodes and edges
    Nnew = insert(i,5);
    Nsimp = Nsimp + Nnew;
    Nnodes = Nnodes + Nnew;
    Nedges = Nedges  + 2*Nnew;
    
    
    % Update edges
    % iedge connects nodes inode1 and Nnodes-Nnew+1
    if(Gamma.edges(iedge,1) == inode1)
        Gamma.edges(iedge,2) = Nnodes - Nnew + 1; % former inode2
    else
        Gamma.edges(iedge,1) = Nnodes - Nnew + 1; % former inode2
    end
    for j=1:Nnew
        Gamma.edges(Nedges - 2*Nnew + 2*(j-1) + 1,:) = [ik0, Nnodes-Nnew+j];
        Gamma.edges(Nedges - 2*Nnew + 2*(j-1) + 2,:) = [Nnodes-Nnew+j, Nnodes-Nnew+j+1]; 
    end
    % Correct last entry
    Gamma.edges(Nedges,:) = [Nnodes,inode2]; 
    
    
    % Compute the new nodes and update Gamma.X
    Xloc_new = zeros(Nnew,3); 
    X1 = Gamma.X(inode1,:); 
    X2 = Gamma.X(inode2,:); 
    
    for j=1:Nnew
        Xloc_new(j,:) = X1 + j/(Nnew+1) * (X2-X1); 
    end
    Gamma.X((Nnodes-Nnew+1):Nnodes,:) = Xloc_new; 
    
    
    
    % Prepare simplex entry for the Nnew new simplices
    for j=1:Nnew
        Gamma.simplices{Nsimp-Nnew+j,1}.nodes = [0,0,0]; 
        Gamma.simplices{Nsimp-Nnew+j,1}.neigh = [0,0,0]; 
        Gamma.simplices{Nsimp-Nnew+j,1}.edges = [0,0,0]; 
        Gamma.simplices{Nsimp-Nnew+j,1}.index = Gamma.simplices{isimp,1}.index; 
    end
    
    
    % Adapt isimp
    Gamma.simplices{isimp,1}.nodes(ind_nr(k0,2)) = Nnodes-Nnew+1; 
    
    neigh_save = Gamma.simplices{isimp,1}.neigh(ind_nr(k0,1));
    Gamma.simplices{isimp,1}.neigh(ind_nr(k0,1)) = Nsimp-Nnew+1;
    
    edge_save = Gamma.simplices{isimp,1}.edges(ind_nr(k0,1)); 
    Gamma.simplices{isimp,1}.edges(ind_nr(k0,1)) = Nedges - 2*Nnew + 1; 
    
    
    % New simplices
    if(Nnew >= 2)
        %j=1;
        Gamma.simplices{Nsimp-Nnew+1,1}.nodes = [Nnodes-Nnew+1,   Nnodes-Nnew+2,   ik0            ];
        Gamma.simplices{Nsimp-Nnew+1,1}.neigh = [Nsimp-Nnew+2,    isimp,           -1             ];
        Gamma.simplices{Nsimp-Nnew+1,1}.edges = [Nedges-2*Nnew+3, Nedges-2*Nnew+1, Nedges-2*Nnew+2];
        %j=2:(Nnew-1)
        for j=2:(Nnew-1)
            Gamma.simplices{Nsimp-Nnew+j,1}.nodes = [Nnodes-Nnew+j,         Nnodes-Nnew+j+1,       ik0                ];
            Gamma.simplices{Nsimp-Nnew+j,1}.neigh = [Nsimp-Nnew+j+1,        Nsimp-Nnew+j-1,        -1                 ];
            Gamma.simplices{Nsimp-Nnew+j,1}.edges = [Nedges-2*Nnew + 2*j+1, Nedges-2*Nnew + 2*j-1, Nedges-2*Nnew + 2*j];
        end
        %j=Nnew
        Gamma.simplices{Nsimp,1}.nodes = [Nnodes,     inode2,   ik0   ];
        Gamma.simplices{Nsimp,1}.neigh = [neigh_save, Nsimp-1,  -1    ];
        Gamma.simplices{Nsimp,1}.edges = [edge_save,  Nedges-1, Nedges];
    else
        if(Nnew == 1)
            Gamma.simplices{Nsimp,1}.nodes = [Nnodes, inode2, ik0];
            Gamma.simplices{Nsimp,1}.neigh = [neigh_save, isimp, -1];  % important here to distinguish case Nnew = 1!
            Gamma.simplices{Nsimp,1}.edges = [edge_save, Nedges-1, Nedges];
        end
    end
    
    % Adapt neigh save
    found = 0; 
    k=1;
    while(k<=3 && found==0)
        if(Gamma.simplices{neigh_save,1}.neigh(k) == isimp)
            found = k;
        end
        k=k+1;
    end
    Gamma.simplices{neigh_save,1}.neigh(found) = Nsimp; 
    
    % Update index, index_edges, index_simp
    index = [index; ((Nnodes-Nnew+1):Nnodes)'];
    index_edges = [index_edges; ((Nedges-2*Nnew+2):2:Nedges)']; 

    index_simp = [index_simp; ((Nsimp-Nnew+1):Nsimp)']; 
    
end



%% Output
Gamma_new = Gamma; 
index_new = index;
index_edges_new = index_edges;
index_simp_new = index_simp; 

% Compute index_edges_corners
M = size(index_edges,1);
index_edges_corners = zeros(M,2);
for k=1:M
    isimp = index_simp(k,1);
    found = 0; 
    j=1; 
    while(j<=3 && found == 0)
        if(Gamma.simplices{isimp,1}.edges(j) == index_edges(k,1))
            found = j;
        end
        j=j+1;
    end
    kloc = found ;
    
    index_edges_corners(k,:) = [Gamma.simplices{isimp,1}.nodes(ind_nr(kloc,1)), Gamma.simplices{isimp,1}.nodes(ind_nr(kloc,2))];

end
end