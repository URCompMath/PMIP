function [X_new,simplices_new,edges_new,index_ext_new, index_nrs] = improve_mesh_after_closing_hole(X,simplices,edges,index_ext, i_p0)
% i_po(r): 1x1, global index of mid point of set
% index_ext: 1x1, one representative for the set

% Initialize output
X_new = X;
simplices_new = simplices;
edges_new = edges;
index_ext_new = index_ext;


nr_simplices = size(simplices,1);

% Indices of new simplices
index_nrs = zeros(1,1);


% Get indices_nrs as well as longest edge information
r=1;

longest_edge_value = 0;
i=index_ext(r);
count = 1;
while(i~= index_ext(r) || count==1)
    % Increase index_nrs(r)
    index_nrs(r) = index_nrs(r)+1;
    
    
    % Get local index of p0
    for i_loc=1:3
        
        if(simplices{i,1}.nodes(i_loc)==i_p0(r))
            i_p0_loc = i_loc;
        end
    end
    
    switch i_p0_loc
        case 1
            ii = 3;
        case 2
            ii = 1;
        case 3
            ii = 2;
    end
    
    % Index of simplex with longest edge (only needed if index_nrs is odd)
    e_global = simplices{i,1}.edges(i_p0_loc);
    length = norm(X(edges(e_global,1),:)-X(edges(e_global,2),:));
    
    if(length > longest_edge_value)
        longest_edge_index_simplex = i;
        longest_edge_value         = length;
        longest_edge_ind_loc       = i_p0_loc;
    end
    
    % Update i, increase count
    count = count + 1;
    i=simplices{i,1}.neigh(ii);
    
end

% Divide simplex with longest edge in 3 new simplices if index_nrs(r) is
% odd

if(index_nrs(r) ~= 2*round(index_nrs(r)/2))  % true --> odd
    nr_simplices = nr_simplices + 2;
    i = longest_edge_index_simplex;
    
    i1 = simplices{i,1}.nodes(1);
    i2 = simplices{i,1}.nodes(2);
    i3 = simplices{i,1}.nodes(3);
    
    % Extend X_new
    p_new = (X(i1,:)+X(i2,:)+X(i3,:))/3;
    X_new = [X_new;p_new];
    X_nr  = size(X_new,1);
    
    
    % local index of p0
    k = longest_edge_ind_loc;
    
    % other local indices
    index = zeros(2,1);
    switch k
        case 1
            index = [3;2];
        case 2
            index = [1;3];
        case 3
            index = [2;1];
    end
    
    % Extend edges_new
    edges_new = [edges_new; zeros(3,2)];
    
    edges_nr = size(edges_new,1);
    edges_new(edges_nr-2,:) = [simplices{i,1}.nodes(index(1)), X_nr];
    edges_new(edges_nr-1,:) = [simplices{i,1}.nodes(index(2)), X_nr];
    edges_new(edges_nr,:)   = [simplices{i,1}.nodes(k), X_nr];
    
    
    % New simplices structure
    simplices_new{i,1} = struct('nodes', [simplices{i,1}.nodes(index(2)), simplices{i,1}.nodes(index(1)), X_nr], ...
        'edges', [edges_nr-2,edges_nr-1,simplices{i,1}.edges(k)],...
        'neigh', [nr_simplices-1, nr_simplices, simplices{i,1}.neigh(k)],...
        'index', simplices{i,1}.index);
    
    sigma1             = struct('nodes', [simplices{i,1}.nodes(index(1)), simplices{i,1}.nodes(k), X_nr], ...
        'edges', [edges_nr,edges_nr-2,simplices{i,1}.edges(index(2))],...
        'neigh', [nr_simplices, i, simplices{i,1}.neigh(index(2))],...
        'index', simplices{i,1}.index);
    
    sigma2             = struct('nodes', [simplices{i,1}.nodes(k), simplices{i,1}.nodes(index(2)), X_nr], ...
        'edges', [edges_nr-1,edges_nr,simplices{i,1}.edges(index(1))],...
        'neigh', [i, nr_simplices-1, simplices{i,1}.neigh(index(1))],...
        'index', simplices{i,1}.index);
    
    % Add to simplices_new
    cell_help = cell(2,1);
    cell_help{1,1} = sigma1;
    cell_help{2,1} = sigma2;
    simplices_new = [simplices_new; cell_help];
    
    % Change neighbor information of neighboring simplices
    i1 = simplices{i,1}.neigh(index(2));
    i2 = simplices{i,1}.neigh(index(1));
    
    e_global_1 = simplices{i,1}.edges(index(2));
    e_global_2 = simplices{i,1}.edges(index(1));
    
    
    for j=1:3
        if(simplices{i1,1}.edges(j)==e_global_1)
            simplices_new{i1,1}.neigh(j) = nr_simplices - 1; % index of sigma 1
        end
        if(simplices{i2,1}.edges(j)==e_global_2)
            simplices_new{i2,1}.neigh(j) = nr_simplices;     % index of simga 2
        end
    end
    
    index_nrs(r) = index_nrs(r) + 1;
    index_ext(r) = nr_simplices;
    
end



% Re-set X, simplices, edges
X = X_new;
simplices = simplices_new;
edges = edges_new;


nr_simplices = size(simplices,1);
nr_edges = size(edges,1);
X_nr     = size(X,1);


% Now, divided 2 neighbor simplices in 4 sub-simplices
r=1;
if(index_nrs(r)>8)
    count = 1;
    i = index_ext(r,1);
    while(count <= index_nrs(r)/2)
        
        % 2 new simplices
        nr_simplices   = nr_simplices + 2;
        
        for j=1:3
            if(simplices{i,1}.nodes(j)==i_p0(r))
                k = j;
            end
        end
        
        index = [0;0];
        switch k
            case 1
                index = [3;2];
            case 2
                index = [1;3];
            case 3
                index = [2;1];
        end
        
        i_plus = simplices{i,1}.neigh(index(1));
        for j=1:3
            if(simplices{i_plus,1}.nodes(j)==i_p0(r))
                k_plus = j;
            end
        end
        index_plus = [0;0];
        switch k_plus
            case 1
                index_plus = [3;2];
            case 2
                index_plus = [1;3];
            case 3
                index_plus = [2;1];
        end
        
        % Calc new point(s) and extend X_new, get index of new nodes, get
        % help neighbor information (new simplices), modify edges_new
        
        if(count == 1)
            p1 = (X(simplices{i,1}.nodes(k),:)+X(simplices{i,1}.nodes(index(1)),:))/2;
            p2 = (X(simplices{i_plus,1}.nodes(k_plus),:)+X(simplices{i_plus,1}.nodes(index_plus(2)),:))/2;
            p = [p1;p2];
            
            X_new = [X_new; p];
            clear p
            
            X_nr  = X_nr + 2;
            
            p_index = [X_nr-1; X_nr];
            p1_help = p_index(1);
            
            aa = nr_simplices + index_nrs(r) -2;
            bb = nr_simplices + 2;
            neigh_help = [aa; bb]; %[nr_simplices + index_nrs(r) -2; nr_simplices + 2];
            
            edges_new(simplices{i,1}.edges(index(2)),:) = [simplices{i,1}.nodes(index(1)), p_index(1)];
            edges_new(simplices{i_plus,1}.edges(index_plus(1)),:) = [simplices{i_plus,1}.nodes(index_plus(2)), p_index(2)];
            
            nr_edges = nr_edges + 4;
            edges_add = [p_index(1), simplices{i,1}.nodes(k);...
                p_index(2), simplices{i,1}.nodes(k);...
                p_index(1), simplices{i,1}.nodes(index(2));...
                p_index(2), simplices{i_plus,1}.nodes(index_plus(1))];
            
            edges_help = [nr_edges-3; nr_edges-2; nr_edges-1; nr_edges];
            edges_down = [nr_edges-3; nr_edges-2];
            
            % Adapt index_ext_new
            index_ext_new(r) = nr_simplices;
            
            
            
        else
            if(count < index_nrs(r)/2)
                p2 = (X(simplices{i_plus,1}.nodes(k_plus),:)+X(simplices{i_plus,1}.nodes(index_plus(2)),:))/2;
                
                
                X_new = [X_new; p2];
                X_nr = X_nr + 1;
                p_index = [X_nr-1; X_nr];
                
                neigh_help = [nr_simplices - 2; nr_simplices + 2];
                
                edges_new(simplices{i_plus,1}.edges(index_plus(1)),:) = [simplices{i_plus,1}.nodes(index_plus(2)), p_index(2)];
                nr_edges = nr_edges + 3;
                edges_add = [p_index(2), simplices{i,1}.nodes(k);...
                    p_index(1), simplices{i,1}.nodes(index(2));...
                    p_index(2), simplices{i_plus,1}.nodes(index_plus(1))];
                
                
                edges_help = [edges_down(2,1); nr_edges-2; nr_edges-1; nr_edges];
                edges_down = [edges_down(1,1); nr_edges-2];
                
            else
                
                p_index = [X_nr; p1_help];
                neigh_help = [nr_simplices - 2; nr_simplices - index_nrs(r) + 2];
                
                nr_edges = nr_edges + 2;
                edges_add = [p_index(1), simplices{i,1}.nodes(index(2));...
                    p_index(2), simplices{i_plus,1}.nodes(index_plus(1))];
                
                edges_help = [edges_down(2,1); edges_down(1,1); nr_edges-1; nr_edges];
            end
        end
        
        
        % Modify edges, calc 3 new edges and extend edges
        
        edges_new(simplices{i,1}.edges(index(1)),:)=[p_index(1), p_index(2)];
        edges_new = [edges_new; edges_add];
        
        
        % New simplices structure
        simplices_new{i,1}      = struct('nodes', [simplices{i,1}.nodes(index(1)), p_index(1), simplices{i,1}.nodes(index(2))], ...
            'edges', [edges_help(3),simplices{i,1}.edges(k),simplices{i,1}.edges(index(2))],...
            'neigh', [nr_simplices-1,simplices{i,1}.neigh(k),simplices{i,1}.neigh(index(2))],...
            'index', simplices{i,1}.index);
        
        simplices_new{i_plus,1} = struct('nodes', [p_index(2), simplices{i_plus,1}.nodes(index_plus(2)), simplices{i_plus,1}.nodes(index_plus(1))], ...
            'edges', [simplices{i_plus,1}.edges(k_plus),edges_help(4),simplices{i_plus,1}.edges(index_plus(1))],...
            'neigh', [simplices{i_plus,1}.neigh(k_plus),nr_simplices-1,simplices{i_plus,1}.neigh(index_plus(1))],...
            'index', simplices{i,1}.index);
        
        sigma1                  = struct('nodes', [p_index(1), p_index(2), simplices{i,1}.nodes(index(2))],...
            'edges', [nr_edges,nr_edges-1,simplices{i,1}.edges(index(1))],...
            'neigh', [i_plus,i,nr_simplices],...
            'index', simplices{i,1}.index);
        
        
        sigma2                  = struct('nodes', [p_index(2), p_index(1), simplices{i,1}.nodes(k)],...
            'edges', [edges_help(1),edges_help(2),simplices{i,1}.edges(index(1))],...
            'neigh', [neigh_help(1),neigh_help(2),nr_simplices-1],...
            'index', simplices{i,1}.index);
        
        
        % Extend simplices_new
        cell_help = cell(2,1);
        cell_help{1,1} = sigma1;
        cell_help{2,1} = sigma2;
        simplices_new  = [simplices_new; cell_help];
        
        % Go to next simplex
        i = simplices{i_plus,1}.neigh(index_plus(1));
        count = count + 1;
        
    end
end



% Divide index_nrs by 2
index_nrs = index_nrs/2;

end