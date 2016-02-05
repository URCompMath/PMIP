function [Gamma_new, top_change] = topological_changes_four_cases(Gamma,config)

Gamma_save = Gamma;
warning = 0; 


%% STEP A: Detect top change Get indices of points where a topology change is likely to occur
size_x = max([abs(min(Gamma.X(:,1))), abs(max(Gamma.X(:,1)))]); 
size_y = max([abs(min(Gamma.X(:,2))), abs(max(Gamma.X(:,2)))]); 
size_z = max([abs(min(Gamma.X(:,3))), abs(max(Gamma.X(:,3)))]); 

sizes = [size_x,size_y,size_z]; 
cc = [-size_x,-size_y,-size_z];

[index_klm,A,Xindex,omega] = get_indices_top_change_four_cases(Gamma,config,config.topchange.a,sizes,cc,config.topchange.Ndetect);


%% STEP B: Identification of the topology change (cases: splitting, merging, increasing genus, decreasing genus)
% Default output (for example if index_klm is empty)
Gamma_new = Gamma;
top_change = 0; 

nr_cubes = size(index_klm,1);
if(nr_cubes == 0)
    return;
end
fprintf('Function: topological changes four cases: proceed\n'); 


% Step B.1 Select cube with largest number of nodes
rank_nr_nodes = zeros(nr_cubes,1); 
for i=1:nr_cubes
    k = index_klm(i,1); 
    l = index_klm(i,2); 
    m = index_klm(i,3); 
    
    rank_nr_nodes(i,1) = A{k,l,m}.nr; 
    
end
[nr_nodes_sort,index_nr_nodes] = sort(rank_nr_nodes); 

if(nr_nodes_sort(nr_cubes) < 2) 
    return;
end
fprintf('Function: topological changes four cases: proceed 2\n'); 


repeat = 0; 
stop = 0; 
while(repeat < nr_cubes && stop == 0)  % idea: if top change == 0, consider next largest cube
    i0 = index_nr_nodes(nr_cubes-repeat,1);
    
    
    % Step B.2 Identification of the top change
    % Create collection of nodes involved in the possible top change
    nodes = zeros(0,1);
    node_nr = 0;
    
    k0 = index_klm(i0,1);
    l0 = index_klm(i0,2);
    m0 = index_klm(i0,3);
    
    
    for j=1:A{k0,l0,m0}.nr
        k_node  = A{k0,l0,m0}.nodes(j,1);
        node_nr = node_nr + 1;
        nodes(node_nr,1) = k_node;
    end
    
    s = 1;
    for alpha=-s:s
        for beta=-s:s
            for gamma=-s:s
                k = min([size(A,1),max([1,k0+alpha])]);
                l = min([size(A,2),max([1,l0+beta])]);
                m = min([size(A,3),max([1,m0+gamma])]);
                
                diff = abs(k-k0)+abs(l-l0)+abs(m-m0); 
                if(diff > 0)
                    % Add (k,l,m) nodes:
                    if(~isempty(A{k,l,m}))
                        for j=1:A{k,l,m}.nr
                            k_node = A{k,l,m}.nodes(j,1);
                            node_nr = node_nr + 1;
                            nodes(node_nr,1) = k_node;
                        end
                    end
                end
            end
        end
    end
                        
    
    % Call fcn identity top change
    top_change_index = identify_top_change(nodes,Xindex,omega,Gamma.X,config);
    
    
    if(top_change_index == 0)
        repeat = repeat + 1;
    else
        stop = 1;
    end
end


fprintf('Top change index %d\n', top_change_index); 
if(top_change_index == 0)
    top_change = 0;
    return;
else
    top_change = 1;
end

%% STEP C: Execution of top change
if(top_change_index == 1 || top_change_index == 4)
        % Splitting or decreasing genus

        % STEP C.1: Determine seperation plane E
        [n0,p0] = get_separation_plane(Gamma.X,nodes);
   
        
        % STEP C.2: Mark simplices for deletion 
        if(top_change_index == 1)
            % Mark simplices for deletion if E \cap sigma \neq
            % \emptyset, get new simplices{i,1}.index (in case of
            % splitting)
            [mark,Gamma.simplices] = mark_simplices_for_deletion(Gamma.X,Gamma.simplices,p0,n0,Gamma.nr_surfaces,config);
            
            Gamma.nr_surfaces = Gamma.nr_surfaces + 1;
        else
            % Decreasing genus: Mark simplices inside a cylinder given by
            % p0 and n0, height 10a and radius 10a
            mark = mark_simplices_for_deletion_cylinder(Gamma.X,Gamma.simplices,p0,n0,config.topchange.a);
        end
          
        
        % STEP C.3: Delete marked simplices, set neigh info of neighbor
        %           simplices to -1
        Gamma = delete_simplices(Gamma,mark);
        fprintf('Plot new surface(s) with temporary holes\n'); 
        figure(1)
        plot_surface(Gamma.X, Gamma.simplices,config.image.flag);
        pause(1)

    
        % STEP C.3b: Mark simplices for deletion of neigh = [-1 -1 -1] or
        % two edges are -1 
        mark = zeros(size(Gamma.simplices,1),1);
        for i=1:size(Gamma.simplices,1)
            neigh = Gamma.simplices{i,1}.neigh;
            Nfree = 0; 
            for j=1:3
                if(neigh(j)==-1)
                    Nfree = Nfree + 1; 
                end
            end
            if(Nfree >= 2)
                mark(i,1)=1;
            end
        end
        
        Gamma = delete_simplices(Gamma,mark);


        
        % STEP C.4: Generate new simplices (close open holes)
        [Gamma.X,Gamma.simplices,Gamma.edges,i_p0] = generate_new_simplices2(Gamma.X,Gamma.simplices,Gamma.edges,p0);
        fprintf('End generate new simplices, before improve mesh\n'); 

        
        % STEP C.5: Improve mesh quality
        for r=1:2
            kmax = 20;
            k=1;
            n_stop = 8; 
            n_remaining = 2*n_stop;  % set it to any value > n_stop
            
            
            while(k<=kmax && n_remaining > n_stop)
                [Gamma,n_remaining] = improve_mesh2(Gamma,i_p0(r));

                k=k+1;
            end
            
        end
                
        Gamma = compute_normal_and_area(Gamma,config);
        figure(1)
        plot_surface(Gamma.X, Gamma.simplices,config.image.flag);
        fprintf('End to splitting/decreasing genus\n'); 
        
        
else
    % Merging or Increasing genus
    if(size(index_klm,1) > 0)

        % STEP C.1: Mark simplices for deletion if part of the collision
        mark = mark_simplices_for_deletion_collision(Gamma.X,Gamma.simplices,nodes);
        Nmark =  sum(mark);

        
        if(Nmark >= 2)
            
            % STEP C.2: Delete marked simplices
            Gamma = delete_simplices(Gamma,mark);
 
            %  Mark simplices for deletion if 2 or 3 edges are free (-1) and
            %  delete them
            mark = zeros(size(Gamma.simplices,1),1);
            for i=1:size(Gamma.simplices,1)
                neigh = Gamma.simplices{i,1}.neigh;
                Nfree = 0;
                for j=1:3
                    if(neigh(j)==-1)
                        Nfree = Nfree + 1;
                    end
                end
                
                if(Nfree >= 2)
                    mark(i,1)=1;
                end
            end
            Gamma = delete_simplices(Gamma,mark);

            
            % Plot intermediate surface with temporary holes
            fprintf('Plot new intermediate surface with temporary holes:\n');
            figure(1)
            plot_surface(Gamma.X, Gamma.simplices, config.image.flag);
            pause(1)
            
            % STEP C.3: Get indices of simplices with one free edge
            [simp_free, v_simp_free, index_simp_free] = get_free_simplices(Gamma);
            % Divide simp_free in two groups, group index stored in column
            % 3
            [simp_free, nr_groups, free_nodes] = get_2_groups(Gamma,simp_free, v_simp_free, index_simp_free);
                            
            if(nr_groups == 2)
                top_change = 1;
                
                % STEP C.4: Get merging pairs
                [Gamma, assign_simp,warning] = get_merging_pairs(Gamma,simp_free);


                % STEP C.5: Merge group1 and group2 according to the assigned
                % pairs, set new surface index for one group
                if(warning == 0)
                    Gamma = merge_edges_after_collision(Gamma, assign_simp);
                end
                fprintf('After merge edges after collision\n');

            else
                if(nr_groups == 1)
                    
                    % STEP C.4: Close open hole, create one new point, connect
                    % all free nodes with it
                    fprintf('One group of free edges (=one hole)\n');
                    Gamma = close_open_hole(Gamma,simp_free,free_nodes);
                    
                    
                    fprintf('Plot new surface\n');
                    plot_surface(Gamma.X, Gamma.simplices, config.image.flag);
                    
                    
                    % STEP C.4b: Improve mesh
                    kmax = 10;
                    k = 1;
                    index_nrs = 100;
                    index_ext = size(Gamma.simplices,1);
                    i_p0 = size(Gamma.X,1);
                    while((index_nrs> 8) && k <= kmax)
                        [Gamma.X,Gamma.simplices,Gamma.edges,index_ext,index_nrs] = improve_mesh_after_closing_hole(Gamma.X,Gamma.simplices,Gamma.edges,index_ext, i_p0);
                        k = k+1;
                                                
                    end
                else
                    fprintf('Error: nr_groups = %d. Cannot be handled yet\n',nr_groups);
                end
            end
            
            
            
        end
        Gamma = compute_normal_and_area(Gamma,config);
        figure(1)
        plot_surface(Gamma.X, Gamma.simplices, config.image.flag);
        fprintf('End of merging / increasing genus\n');

    end
        
end
        
fprintf('End of top change main fcn\n'); 

 
%% Set output
Gamma_new = Gamma;
if(warning)
    Gamma_new = Gamma_save; 
    top_change = 0; 
    fprintf('Warning active, reset old surfaces!\n'); 
end

end

 