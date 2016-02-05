function [cp,i_simp_closest] = project2surface_parent(p,v,i_parent,Surface)
% fcn close to "get_closest_simp_use_parents", here: is_inside is replaced
% by is_inside_with_dir since we project not orthogonally, but in direction
% of a given direction v

% Consider first parent simplex
[inside,q,delta_lambda] = is_inside_with_dir(p,v,i_parent,Surface); 

if(inside)
    cp = q;
    i_simp_closest = i_parent;
else
    % Save delta_lambda and projected node
    delta_lambda_min = delta_lambda;
    
    % Consider neighbor triangles
    Nsimp = size(Surface.tri,1);
    
    used = zeros(Nsimp,1);
    found = 0;
    found_q = zeros(3,1);
    
    level = 1;
    
    list = Surface.neigh(i_parent,:)';
    
    while(found==0 && level <= Surface.level_max)
        % Set used information
        used(list,1) = 1;
        
        
        % Search in list
        Nlist = size(list,1);
        k=1;
        while(k<= Nlist && found==0)
            i_simp_list = list(k,1);
            [inside,q,delta_lambda] = is_inside_with_dir(p,v,i_simp_list,Surface);
            
            if(inside)
                found = i_simp_list;
                found_q = q;
            else
                if(delta_lambda < delta_lambda_min)
                    delta_lambda_min = delta_lambda;
                end
                k=k+1;
            end
        end
        
        % Consider neighbor from list members if found==0
        if(found == 0)
            list_save = list;
            Nnew = 0;
            list = zeros(Nlist,1);
            for k=1:Nlist
                for l=1:3
                    i_neigh = Surface.neigh(list_save(k,1),l);
                    if(used(i_neigh,1) == 0)
                        Nnew = Nnew + 1;
                        list(Nnew,1) = i_neigh;
                    end
                end
            end
            list = list(1:Nnew,1);
        end
        
        level = level + 1;
        
    end
    if(found > 0)
        i_simp_closest = found;
        cp = found_q;
    else
        fprintf('fcn: project2surface_parent:\nCannot find a closest simplex for p = (%3.5f, %3.5f, %3.5f) near i_parent = %d\n', p(1,1),p(1,2),p(1,3),i_parent);
        pause
    end
end

end

