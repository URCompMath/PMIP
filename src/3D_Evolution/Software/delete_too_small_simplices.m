function [X_new,simplices_new,edges_new] = delete_too_small_simplices(X,simplices,edges,mark,involved_simplices)

% Initialize output
simplices_new = simplices;
X_new = X;
edges_new = edges;

% Sizes
Nsimp = size(simplices,1);
Nnodes = size(X,1);
Nedges = size(edges,1);

% Information about simplices, nodes and edges too be deleted
info_simplices = zeros(Nsimp,1);
info_nodes     = zeros(Nnodes,1);
info_edges     = zeros(Nedges,1);


fprintf('%d simplices marked for possible deletion\n', size(mark,1));
Ndelete = 0;

for count = 1:size(mark,1);
    i = mark(count,1);
    if(i>0)
        delete_flag = mark(count,2);
        if(delete_flag == 1)  % too small area
            % Get nodes which shrink to one node
            j1 = simplices{i,1}.nodes(1);
            j2 = simplices{i,1}.nodes(2);
            j3 = simplices{i,1}.nodes(3);
            
            
            % Get neighbor simplices which are also deleted
            n1 = simplices{i,1}.neigh(1);
            n2 = simplices{i,1}.neigh(2);
            n3 = simplices{i,1}.neigh(3);
            
            % Stop if i, n1, n2 or n3 is a simplex involved in previous
            % refinement
            stop = 0;
            if(involved_simplices(i,1)==1 || involved_simplices(n1,1)==1 || involved_simplices(n2,1)==1 || involved_simplices(n3,1)==1)
                stop = 1;
            end
            
            if(stop == 0)
                Ndelete = Ndelete + 4;
                
                % Get edges which to be deleted
                e1 = simplices{i,1}.edges(1);
                e2 = simplices{i,1}.edges(2);
                e3 = simplices{i,1}.edges(3);
                
                % Get neighbors and edges of neighbors
                j11_loc = find(simplices{n1,1}.neigh == i);
                j12_loc = mod(j11_loc,3)+1;
                j13_loc = mod(j12_loc,3)+1;
                n11 = simplices{n1,1}.neigh(j12_loc);
                n12 = simplices{n1,1}.neigh(j13_loc);
                e11 = simplices{n1,1}.edges(j12_loc);
                e12 = simplices{n1,1}.edges(j13_loc);
                
                j21_loc = find(simplices{n2,1}.neigh == i);
                j22_loc = mod(j21_loc,3)+1;
                j23_loc = mod(j22_loc,3)+1;
                n21 = simplices{n2,1}.neigh(j22_loc);
                n22 = simplices{n2,1}.neigh(j23_loc);
                e21 = simplices{n2,1}.edges(j22_loc);
                e22 = simplices{n2,1}.edges(j23_loc);
                
                j31_loc = find(simplices{n3,1}.neigh == i);
                j32_loc = mod(j31_loc,3)+1;
                j33_loc = mod(j32_loc,3)+1;
                n31 = simplices{n3,1}.neigh(j32_loc);
                n32 = simplices{n3,1}.neigh(j33_loc);
                e31 = simplices{n3,1}.edges(j32_loc);
                e32 = simplices{n3,1}.edges(j33_loc);
                
                % If nj1 and nj2 are already neighbors, i is not deleted
                found0 = 0;
                help_neigh = simplices{n11,1}.neigh;
                if(~isempty(find(help_neigh == n12)))
                    found0 = 1;
                end
                help_neigh = simplices{n21,1}.neigh;
                if(~isempty(find(help_neigh == n22)))
                    found0 = 1;
                end
                help_neigh = simplices{n31,1}.neigh;
                if(~isempty(find(help_neigh == n32)))
                    found0 = 1;
                end
                
                % Proceed if found == 0
                if(found0 == 0)
                    
                    special_case = 0;
                    if(n11==n3 || n12==n2 || n21==n1 || n22==n3 || n31==n2 || n32==n1)
                        special_case = 1;
                    end
                    
                    if(special_case == 0)
                        % Set flag for deletion of nodes, simplices and edges
                        info_simplices(i,1) = 1;
                        info_simplices(n1,1) = 1;
                        info_simplices(n2,1) = 1;
                        info_simplices(n3,1) = 1;
                        
                        info_nodes(j2,1) = 1;
                        info_nodes(j3,1) = 1;
                        
                        info_edges(e1,1)  = 1;
                        info_edges(e2,1)  = 1;
                        info_edges(e3,1)  = 1;
                        info_edges(e12,1) = 1;
                        info_edges(e22,1) = 1;
                        info_edges(e32,1) = 1;
                        
                        % Adaptations
                        % Shrink j1, j2 and j3 to one node
                        X_new(j1,:) = (X(j1,:) + X(j2,:) + X(j3,:))/3;
                        
                        % Change node info of all simplices containing vertex j2 or j3
                        % to j1, set mark(:,1) to 0 if simplex is in list for deletion
                        for isimp = 1:Nsimp
                            nodes_help = simplices{isimp,1}.nodes;
                            k = find(nodes_help == j2);
                            if(isempty(k)==0)
                                simplices_new{isimp,1}.nodes(k) = j1;
                                l = find(mark(:,1) == isimp);
                                if(isempty(l)==0)
                                    mark(l,1) = 0;
                                end
                            end
                            k = find(nodes_help == j3);
                            if(isempty(k)==0)
                                simplices_new{isimp,1}.nodes(k) = j1;
                                l = find(mark(:,1) == isimp);
                                if(isempty(l)==0)
                                    mark(l,1) = 0;
                                end
                            end
                            k = find(nodes_help == j1);
                            if(isempty(k)==0)
                                l = find(mark(:,1) == isimp);
                                if(isempty(l)==0)
                                    mark(l,1) = 0;
                                end
                            end
                            
                            
                            
                        end
                        % Set mark(:,1) to 0 for neighbors
                        n_loop = 8;
                        index_loop = simplices_new{isimp,1}.neigh';
                        index_loop_new = zeros(0,1);
                        n_new = 0;
                        for i_loop = 1:n_loop
                            for ii = 1: size(index_loop,1)
                                l = find(mark(:,1) == index_loop(ii,1));
                                if(isempty(l)==0)
                                    mark(l,1) = 0;
                                end
                                n_new = n_new + 3;
                                index_loop_new((n_new-2):n_new,1) = simplices{index_loop(ii,1),1}.neigh';
                            end
                            index_loop = index_loop_new;
                            index_loop_new = zeros(0,1);
                            n_new = 0;
                        end
                        
                        
                        % Change node info of all edges containing nodes j2 or j3 to j1
                        for iedges = 1:Nedges
                            k = find(edges(iedges,:) == j2);
                            if(isempty(k)==0)
                                edges_new(iedges,k) = j1;
                            else
                                k2 = find(edges(iedges,:) == j3);
                                if(isempty(k2)==0)
                                    edges_new(iedges,k2) = j1;
                                end
                            end
                        end
                        
                        % Adapt edge and neigh info of neighbors n11-n31, nodes already adapted
                        % n11 (edges ok, change neighbor from n1 to n12)
                        neigh_help = simplices{n11,1}.neigh;
                        simplices_new{n11,1}.neigh(find(neigh_help == n1)) = n12;
                        
                        % n12 (change edge from e12 to e11, neighbor from n1 to n11)
                        edges_help = simplices{n12,1}.edges;
                        simplices_new{n12,1}.edges(find(edges_help == e12)) = e11;
                        neigh_help = simplices{n12,1}.neigh;
                        simplices_new{n12,1}.neigh(find(neigh_help == n1)) = n11;
                        
                        % n21 (edges ok, change neighbor form n2 to n22)
                        neigh_help = simplices{n21,1}.neigh;
                        simplices_new{n21,1}.neigh(find(neigh_help == n2)) = n22;
                        
                        % n22 (change edge from e22 to e21, neighbor from n2 to n21)
                        edges_help = simplices{n22,1}.edges;
                        simplices_new{n22,1}.edges(find(edges_help == e22)) = e21;
                        neigh_help = simplices{n22,1}.neigh;
                        simplices_new{n22,1}.neigh(find(neigh_help == n2)) = n21;
                        
                        % n31 (edges ok, change neighbor form n3 to n32)
                        neigh_help = simplices{n31,1}.neigh;
                        simplices_new{n31,1}.neigh(find(neigh_help == n3)) = n32;
                        
                        % n32 (change edge from e32 to e31, neighbor from n3 to n31)
                        edges_help = simplices{n32,1}.edges;
                        simplices_new{n32,1}.edges(find(edges_help == e32)) = e31;
                        neigh_help = simplices{n32,1}.neigh;
                        simplices_new{n32,1}.neigh(find(neigh_help == n3)) = n31;
                        
                        
                    else
                        % Possible shift of indices
                        if(n22==n3)  % j2 -> j1, j3 -> j2
                            help = j1;
                            j1 = j2;
                            j2 = j3;
                            j3 = help;
                            
                            help = n1;
                            n1 = n2;
                            n2 = n3;
                            n3 = help;
                            
                            help = e1;
                            e1 = e2;
                            e2 = e3;
                            e3 = help;
                                                        
                            help = n11;
                            n11 = n21;
                            n21 = n31;
                            n31 = help;
                            
                            help = n12;
                            n12 = n22;
                            n22 = n32;
                            n32 = help;
                            
                            help = e11;
                            e11 = e21;
                            e21 = e31;
                            e31 = help;
                            
                            help = e12;
                            e12 = e22;
                            e22 = e32;
                            e32 = help;
                            
                        else
                            if(n32==n1)% j3 -> j1, j2 -> j3
                                help = j1;
                                j1 = j3;
                                j3 = j2;
                                j2 = help;
                                
                                help = n1;
                                n1 = n3;
                                n3 = n2;
                                n2 = help;
                                
                                help = e1;
                                e1 = e3;
                                e3 = e2;
                                e2 = help;
                                                                
                                help = n11;
                                n11 = n31;
                                n31 = n21;
                                n21 = help;
                                
                                help = n12;
                                n12 = n32;
                                n32 = n22;
                                n22 = help;
                                
                                help = e11;
                                e11 = e31;
                                e31 = e21;
                                e21 = help;
                                
                                help = e12;
                                e12 = e32;
                                e32 = e22;
                                e22 = help;
                                
                            end
                        end
                        % Now we have the situation in which n12 == n2
                        % and n21 == n1
                        
                        % Do not delete in special cases
                        found1 = 0;
                        if(n22==n3 || n11==n3)
                            found1 = 1;
                        end
                        
                        % If one of the nij are already neighbors, do not
                        % delete i
                        help_neigh = simplices{n11,1}.neigh;
                        if(~isempty(find(help_neigh == n22)))
                            found1 = 1;
                        end
                        help_neigh = simplices{n31,1}.neigh;
                        if(~isempty(find(help_neigh == n32)))
                            found1 = 1;
                        end
                        
                        % Proceed only if found1 == 0
                        if(found1 == 0)
                            
                            % Set flag for deletion of nodes, simplices and edges
                            info_simplices(i,1) = 1;
                            info_simplices(n1,1) = 1;
                            info_simplices(n2,1) = 1;
                            info_simplices(n3,1) = 1;
                            
                            info_nodes(j2,1) = 1;
                            info_nodes(j3,1) = 1;
                            
                            info_edges(e1,1)  = 1;
                            info_edges(e2,1)  = 1;
                            info_edges(e3,1)  = 1;
                            info_edges(e12,1) = 1; % note e12 == e21!!!, e21 is also deleted here!!!, e11 and e31 remain
                            info_edges(e22,1) = 1;
                            info_edges(e32,1) = 1;
                            
                            % Adaptations
                            % Shrink j1, j2 and j3 to one node
                            X_new(j1,:) = (X(j1,:) + X(j2,:))/2; % use j1 and j2 only!
                            
                            % Change node info of all simplices containing vertex j2 or j3
                            % to j1, set mark(:,1) to 0 if simplex is in list for deletion
                            for isimp = 1:Nsimp
                                nodes_help = simplices{isimp,1}.nodes;
                                k = find(nodes_help == j2);
                                if(isempty(k)==0)
                                    simplices_new{isimp,1}.nodes(k) = j1;
                                    l = find(mark(:,1) == isimp);
                                    if(isempty(l)==0)
                                        mark(l,1) = 0;
                                    end
                                end
                                k = find(nodes_help == j3);
                                if(isempty(k)==0)
                                    simplices_new{isimp,1}.nodes(k) = j1;
                                    l = find(mark(:,1) == isimp);
                                    if(isempty(l)==0)
                                        mark(l,1) = 0;
                                    end
                                end
                                k = find(nodes_help == j1);
                                if(isempty(k)==0)
                                    l = find(mark(:,1) == isimp);
                                    if(isempty(l)==0)
                                        mark(l,1) = 0;
                                    end
                                end
                            end
                            % Change node info of all edges containing nodes j2 or j3 to j1
                            for iedges = 1:Nedges
                                k = find(edges(iedges,:) == j2);
                                if(isempty(k)==0)
                                    edges_new(iedges,k) = j1;
                                else
                                    k2 = find(edges(iedges,:) == j3);
                                    if(isempty(k2)==0)
                                        edges_new(iedges,k2) = j1;
                                    end
                                end
                            end

                            
                            % Adapt edge and neigh info of neighbors n11-n31, nodes already adapted
                            % n11 (edges ok, change neighbor from n1 to n22)
                            % Attention: n12 == n2!
                            neigh_help = simplices{n11,1}.neigh;
                            simplices_new{n11,1}.neigh(find(neigh_help == n1)) = n22;
                            
                            % n22 (change edge from e22 to e11, neighbor from n2 to n11)
                            % Attention: n21 == n1!
                            edges_help = simplices{n22,1}.edges;
                            simplices_new{n22,1}.edges(find(edges_help == e22)) = e11;
                            neigh_help = simplices{n22,1}.neigh;
                            simplices_new{n22,1}.neigh(find(neigh_help == n2)) = n11;
                            
                            % n31 (edges ok, change neighbor form n3 to n32)
                            neigh_help = simplices{n31,1}.neigh;
                            simplices_new{n31,1}.neigh(find(neigh_help == n3)) = n32;
                            
                            % n32 (change edge from e32 to e31, neighbor from n3 to n31)
                            edges_help = simplices{n32,1}.edges;
                            simplices_new{n32,1}.edges(find(edges_help == e32)) = e31;
                            neigh_help = simplices{n32,1}.neigh;
                            simplices_new{n32,1}.neigh(find(neigh_help == n3)) = n31;
                            

                        end
                    end
                    
                    
                    
                    
                    % Update simplices, X, edges
                    simplices = simplices_new;
                    edges = edges_new;
                    X = X_new;
                end
            end
        else
            % One small edge, two larger edges
            % Find small edge
            j1_temp = simplices{i,1}.nodes(1);
            j2_temp = simplices{i,1}.nodes(2);
            j3_temp = simplices{i,1}.nodes(3);
            
            X1 = X(j1_temp,:);
            X2 = X(j2_temp,:);
            X3 = X(j3_temp,:);
            
            h1 = norm(X3-X2);
            h2 = norm(X1-X3);
            h3 = norm(X2-X1);
            [h,index] = min([h1,h2,h3]);
            
            % Local and global node indices
            j1_loc = index;
            j2_loc = mod(j1_loc,3)+1;
            j3_loc = mod(j2_loc,3)+1;
            
            j1 = simplices{i,1}.nodes(j1_loc);
            j2 = simplices{i,1}.nodes(j2_loc);
            j3 = simplices{i,1}.nodes(j3_loc); % j2 and j3 shrink to one node
            
            
            % Get neighbor simplex which is also deleted
            in = simplices{i,1}.neigh(j1_loc);
            
            % Stop if i or in are involved in previous refinement
            % refinement
            stop = 0;
            if(involved_simplices(i,1)==1 || involved_simplices(in,1)==1)
                stop = 1;
            end
            
            if(stop == 0)
                Ndelete = Ndelete + 2;
                
                neigh_help = simplices{in,1}.neigh;
                jn1_loc = find(neigh_help == i);
                jn2_loc = mod(jn1_loc,3)+1;
                jn3_loc = mod(jn2_loc,3)+1;
                
                j4 = simplices{in,1}.nodes(jn1_loc);
                
                % Get edges
                e1 = simplices{i,1}.edges(j1_loc);
                e2 = simplices{i,1}.edges(j2_loc);
                e3 = simplices{i,1}.edges(j3_loc);
                e4 = simplices{in,1}.edges(jn2_loc);
                e5 = simplices{in,1}.edges(jn3_loc);
                
                
                % Get neighbors
                n1 = simplices{i,1}.neigh(j2_loc);
                n2 = simplices{i,1}.neigh(j3_loc);
                n3 = simplices{in,1}.neigh(jn2_loc);
                n4 = simplices{in,1}.neigh(jn3_loc);
                
                % If n1 and n2 or n3 and n4 are already neighbors, do not
                % delete i
                found0 = 0;
                help_neigh = simplices{n1,1}.neigh;
                if(~isempty(find(help_neigh == n2)))
                    found0 = 1;
                end
                help_neigh = simplices{n3,1}.neigh;
                if(~isempty(find(help_neigh == n4)))
                    found0 = 1;
                end
                
                % Proceed only if found0==0
                if(found0 == 0)
                    
                    % Set flag for deletion of nodes, simplices and edges
                    info_simplices(i,1) = 1;
                    info_simplices(in,1) = 1;
                    
                    info_nodes(j3,1) = 1;
                    
                    info_edges(e1,1)  = 1;
                    info_edges(e3,1)  = 1;
                    info_edges(e5,1)  = 1;
                    
                    
                    % Adaptations
                    % Shrink nodes j2 and j3
                    X_new(j2,:) = (X(j2,:) + X(j3,:))/2;
                    
                    % Change node info of all simplices containing vertex j3
                    % to j2, set mark(:,1) to 0 if simplex is in list for deletion
                    for isimp = 1:Nsimp
                        nodes_help = simplices{isimp,1}.nodes;
                        k = find(nodes_help == j3);
                        if(isempty(k)==0)
                            simplices_new{isimp,1}.nodes(k) = j2;
                            l = find(mark(:,1) == isimp);
                            if(isempty(l)==0)
                                mark(l,1) = 0;
                            end
                        end
                        k = find(nodes_help == j1);
                        if(isempty(k)==0)
                            l = find(mark(:,1) == isimp);
                            if(isempty(l)==0)
                                mark(l,1) = 0;
                            end
                        end
                        k = find(nodes_help == j2);
                        if(isempty(k)==0)
                            l = find(mark(:,1) == isimp);
                            if(isempty(l)==0)
                                mark(l,1) = 0;
                            end
                        end
                        k = find(nodes_help == j4);
                        if(isempty(k)==0)
                            l = find(mark(:,1) == isimp);
                            if(isempty(l)==0)
                                mark(l,1) = 0;
                            end
                        end
                        
                    end
                    
                    
                    % Change node info of all edges containing nodes j3 to j2
                    for iedges = 1:Nedges
                        k = find(edges(iedges,:) == j3);
                        if(isempty(k)==0)
                            edges_new(iedges,k) = j2;
                        end
                    end
                    
                    % Adapt edge and neigh info of neighbors n1-n4, nodes already adapted
                    % n1 (edges ok, change neighbor from i to n2
                    neigh_help = simplices{n1,1}.neigh;
                    simplices_new{n1,1}.neigh(find(neigh_help == i)) = n2;
                    
                    % n2 (change edge from e3 to e2, neighbor from i to n1)
                    edges_help = simplices{n2,1}.edges;
                    simplices_new{n2,1}.edges(find(edges_help == e3)) = e2;
                    neigh_help = simplices{n2,1}.neigh;
                    simplices_new{n2,1}.neigh(find(neigh_help == i)) = n1;
                    
                    % n3 (edges ok, change neighbor from in to n4
                    neigh_help = simplices{n3,1}.neigh;
                    simplices_new{n3,1}.neigh(find(neigh_help == in)) = n4;
                    
                    % n4 (change edge from e5 to e4, neighbor from in to n3)
                    edges_help = simplices{n4,1}.edges;
                    simplices_new{n4,1}.edges(find(edges_help == e5)) = e4;
                    neigh_help = simplices{n4,1}.neigh;
                    simplices_new{n4,1}.neigh(find(neigh_help == in)) = n3;
                    
                    % Update simplices, X, edges
                    simplices = simplices_new;
                    edges = edges_new;
                    X = X_new;
                end
            end
        end
        
        
    end
end

% Delete simplices, nodes, edges according to info_simplices, info_nodes,
% info_edges and shift other indices

[nodes_del,help]  = find(info_nodes==1);
[sim_del,help]    = find(info_simplices==1);
[edges_del,help]  = find(info_edges==1);



% Shift remaining indices in simplices, nodes, neigh, edges
for i=1:Nsimp
    if(info_simplices(i,1)== 0)
        nodes0 = simplices_new{i,1}.nodes;
        neigh0 = simplices_new{i,1}.neigh;
        edges0 = simplices_new{i,1}.edges;
        
        for j=1:3
            simplices_new{i,1}.nodes(j) = nodes0(j) - size(find(nodes_del < nodes0(j)),1);
            simplices_new{i,1}.neigh(j) = neigh0(j) - size(find(sim_del   < neigh0(j)),1);
            simplices_new{i,1}.edges(j) = edges0(j) - size(find(edges_del < edges0(j)),1);
        end
    end
end
for i=1:Nedges
    if(isempty(find(edges_del == i)))
        for j=1:2
            edges_new(i,j) = edges(i,j) - size(find(nodes_del < edges(i,j)),1);
        end
    end
end
fprintf('Deleted %d simplices (incl. neighbors)\n\n', Ndelete);


% Delete rows in simplices_new and X_new
simplices_new(sim_del,:)=[];
X_new(nodes_del,:)=[];
edges_new(edges_del,:)=[];


end

