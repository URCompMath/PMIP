function [X_new,simplices_new,edges_new,involved_simplices] = refine_marked_simplices(X,simplices,edges,mark)
% mark: N_mark x 1 matrix, contains indexes of simplices to be refined
N_mark = size(mark,1);
fprintf('%d simplices marked for possible refinement, %d refinements possible (incl. neighbors)\n', N_mark, 2*N_mark); 

% Initialize output
simplices_new = simplices;
size_X = size(X,1);
size_edges = size(edges,1);


% Involved Simplices, store involved (refined simplices) such that they are
% not coarsened in the next step
involved_simplices = zeros(2*size(simplices,1),1); 

Nrefine = 0; 
for count=1:N_mark
     
    
    i = mark(count);
    if(i>0)        
        % Search for biggest angle
        v1 = X(simplices{i,1}.nodes(1),:)-X(simplices{i,1}.nodes(2),:);
        v2 = X(simplices{i,1}.nodes(2),:)-X(simplices{i,1}.nodes(3),:);
        v3 = X(simplices{i,1}.nodes(3),:)-X(simplices{i,1}.nodes(1),:);
        a1 = acos(-v1/norm(v1) * (v3/norm(v3))');
        a2 = acos(-v2/norm(v2) * (v1/norm(v1))');
        a3 = acos(-v3/norm(v3) * (v2/norm(v2))');
            
        [a,i3] = max([a1,a2,a3]);
        
        i_refinement_edge = simplices{i,1}.edges(i3);
        i_neigh = simplices{i,1}.neigh(i3);
        [i1,i2] = get_ind(i3);

        % If i_neigh is in list mark, set corresponding mark index to 0
        k = count+1;
        found = 0;
        while(k<=N_mark && found==0)
            if(mark(k)==i_neigh)
                mark(k)=0;
                found = 1;
            end
            k=k+1;
        end
        
        
        % Search for edge index of i_neigh such that i_refinement_edge is
        % the common edge
        found = 0;
        count0 = 1;
        while(count0<=3 && found==0)
            if(i_refinement_edge == simplices{i_neigh,1}.edges(count0))
                found = count0;
            end
            count0 = count0+1;
        end
        j3 = found;
        [j1,j2] = get_ind(j3);
        
            
        % Add 2 new simplexes to simplices_new
        N_new = size(simplices_new,1) + 2 ;
        simplices_new{N_new-1,1} = struct('nodes', [0 0 0], 'edges', [0 0 0], 'neigh', [0 0 0], 'index', simplices{i,1}.index);
        simplices_new{N_new,1}   = struct('nodes', [0 0 0], 'edges', [0 0 0], 'neigh', [0 0 0], 'index', simplices{i,1}.index);
        
        % Create new node on refinement edge
        p = 0.5*(X(simplices{i,1}.nodes(i1),:)+X(simplices{i,1}.nodes(i2),:));
        
        % Add new node to X and update mark_edges
        size_X = size_X + 1;
        i_new  = size_X;
        X(i_new,:) = p;
        
        
        % Create 3 new edges and change edge(i_refinement_edge,:)
        size_edges = size_edges + 3;
        edges(size_edges-2,:)      = [i_new,simplices{i,1}.nodes(i3)];
        edges(size_edges-1,:)      = [i_new,simplices{i,1}.nodes(i2)];
        edges(size_edges,:)        = [i_new,simplices{i_neigh,1}.nodes(j3)];
        edges(i_refinement_edge,:) = [i_new,simplices{i,1}.nodes(i1)];
        
        
        % Set simplices{*,1}.edges
        simplices_new{i,1}.edges       = [i_refinement_edge,size_edges-2,simplices{i,1}.edges(i2)];
        simplices_new{N_new-1,1}.edges = [size_edges-2,size_edges-1,simplices{i,1}.edges(i1)];
        simplices_new{i_neigh,1}.edges = [size_edges,i_refinement_edge,simplices{i_neigh,1}.edges(j1)];
        simplices_new{N_new,1}.edges   = [size_edges-1,size_edges,simplices{i_neigh,1}.edges(j2)];
        
        
        
        % Set simplices_new{*,1}.nodes
        simplices_new{i,1}.nodes       = [simplices{i,1}.nodes(i3), simplices{i,1}.nodes(i1),i_new];
        simplices_new{N_new-1,1}.nodes = [simplices{i,1}.nodes(i2), simplices{i,1}.nodes(i3),i_new];
        simplices_new{i_neigh,1}.nodes = [simplices{i_neigh,1}.nodes(j2), simplices{i_neigh,1}.nodes(j3),i_new];
        simplices_new{N_new,1}.nodes   = [simplices{i_neigh,1}.nodes(j3), simplices{i_neigh,1}.nodes(j1),i_new];
        
        
        % Set simplices_new{i,1}.neigh and simplices_new{N_new,1}.neigh
        simplices_new{i,1}.neigh       = [i_neigh, N_new-1, simplices{i,1}.neigh(i2)];
        simplices_new{N_new-1,1}.neigh = [i, N_new, simplices{i,1}.neigh(i1)];
        simplices_new{i_neigh,1}.neigh = [N_new, i, simplices{i_neigh,1}.neigh(j1)];
        simplices_new{N_new,1}.neigh   = [N_new-1, i_neigh, simplices{i_neigh,1}.neigh(j2)];
        
        % Change neighbor information of neighbors
        for k=1:2
            i_n = simplices_new{N_new-2+k,1}.neigh(3);
            if(i_n>0)
                edge_index = simplices_new{N_new-2+k,1}.edges(3);
                found = 0;
                count0 = 1;
                while(count0<=3 && found==0)
                    if(edge_index == simplices{i_n,1}.edges(count0))
                        found = count0;
                    end
                    count0 = count0+1;
                end
                simplices_new{i_n,1}.neigh(found)= N_new-2+k;
            end
        end
        
        % Update simplices
        simplices = simplices_new;
        
        % Update involved_simplices
        involved_simplices(i,1)       = 1;
        involved_simplices(i_neigh,1) = 1;
        involved_simplices(N_new-1,1) = 1;
        involved_simplices(N_new,1)   = 1;
        
        
        Nrefine = Nrefine + 2; 

    end
end
fprintf('Refined %d simplices (incl. neighbor)\n', Nrefine); 

% Set output
edges_new = edges;
X_new = X;
involved_simplices = involved_simplices(1:size(simplices_new,1),1); 

end


function [j1,j2] = get_ind(j)
switch j
    case 1
        j1 = 2;
        j2 = 3;
    case 2
        j1 = 3;
        j2 = 1;
    case 3
        j1 = 1;
        j2 = 2;
end
end