function [simp_free_new, nr_groups, free_nodes] = get_2_groups(Gamma,simp_free, v_simp_free, index_simp_free)

N = size(simp_free,1); 

% Get group1 and group2

% First option: look for two different surface indices - decommented:
% search is easier, but degenerated cases (e.g. 3 holes) are not detected! 
index1 = index_simp_free(1,1);
simp_free(1,3) = 1; 

index2 = -1;



index2 = -1; 

ind_nr = [2,3; 3,1; 1,2]; 

free_edges = zeros(N,1);
free_nodes = zeros(N,2);

if(index2 == -1)
    % Reset simp_free(:,3)
    simp_free(:,3) = zeros(N,1); 
        
    for i=1:N
        i0 = simp_free(i,1); 
        k0 = simp_free(i,2); 
        free_edges(i,1) = Gamma.simplices{i0,1}.edges(k0); 
        
        nodes = Gamma.simplices{i0,1}.nodes;
        free_nodes(i,:) = [nodes(ind_nr(k0,1)), nodes(ind_nr(k0,2))]; 
    end
    
    % 1st loop
    istart = free_nodes(1,1);
    simp_free(1,3) = 1; 
    
    i = free_nodes(1,2); 
    search = free_nodes(:,1); 
    while(i ~= istart)
        k = find(search == i); 
        if(isempty(k))
            fprintf('Error in get_2_groups, probably orientation error in free nodes or corrupted mesh\n'); 
            pause
        else
            i = free_nodes(k,2); 
            simp_free(k,3) = 1; 
        end
    end
    
    % 2nd loop if k_group2 is not empty
    search = simp_free(:,3); 
    k_group2 = find(search == 0);
    if(isempty(k_group2))
        fprintf('One group\n'); 
        nr_groups = 1; 
    else
        fprintf('Two groups\n'); 
        nr_groups = 2; 
        
        istart = free_nodes(k_group2(1),1); 
        simp_free(k_group2(1),3) = 2; 
        
        i = free_nodes(k_group2(1),2); 
        search = free_nodes(:,1); 
        
        while(i~= istart)
            k = find(search == i);
            if(isempty(k))
                fprintf('Error in get_2_groups, probably orientation error in free nodes or corrupted mesh\n');
                pause
            else
                i = free_nodes(k,2); 
                simp_free(k,3) = 2;
            end
        end
    end
       
    % Test for 3rd group
    search = simp_free(:,3); 
    k_group2 = find(search == 0);
    if(~isempty(k_group2))
        fprintf('Error: At least three open holes occur. Cannot be handled yet\n'); 
    end
    
    

end

% Get output
simp_free_new = simp_free;
end
    
