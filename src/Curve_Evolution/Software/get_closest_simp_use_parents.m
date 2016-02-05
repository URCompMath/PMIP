function [i_simp_closest,cp] = get_closest_simp_use_parents(p,i_parent,Surface)

% Consider first parent simplex
[inside,q,delta_lambda] = is_inside(p,i_parent,Surface);

if(inside)
    cp = q;
    i_simp_closest = i_parent;
else
    % Save delta_lambda and projected node
    delta_lambda_min = delta_lambda;
    q_min = q;
    i_min = i_parent;
    
    % Look for neighbors
    Nsimp = size(Surface.tri,1);
    
    used = zeros(Nsimp,1);
    found = 0;
    found_q = zeros(3,1);
    
    level = 1;
    
    list = Surface.neigh(i_parent,:)';
    
    while(found==0 && level <= Surface.level_max)
        % Set used information
        used(list,1) = 1;
        
        % Set improve flag
        improve_flag = 0;

        % Search in list
        Nlist = size(list,1);
        k=1;
        while(k<= Nlist && found==0)
            i_simp_list = list(k,1);
            [inside,q] = is_inside(p,i_simp_list,Surface);
            
            if(inside)
                found = i_simp_list;
                found_q = q;
            else
                if(delta_lambda < delta_lambda_min)
                    delta_lambda_min = delta_lambda;
                    q_min = q;
                    i_min = i_simp_list;
                    improve_flag = 1;
                end
                k=k+1;
            end
        end
        
        % Consider neighbor from list members if found==0
        if(found == 0)
            if(improve_flag == 1)
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
            else
                % No improvement, use former projection and simplex
                found = i_min; 
                found_q = q_min;
            end
        end
        
        level = level + 1;
        
    end
    if(found > 0)
        i_simp_closest = found;
        cp = found_q;
    else
        fprintf('fcn: get_closest_simp_use_parents:\nCannot find a closest simplex for p = (%3.5f, %3.5f, %3.5f) near i_parent = %d\n', p(1,1),p(1,2),p(1,3),i_parent);
        pause
    end
end

end
