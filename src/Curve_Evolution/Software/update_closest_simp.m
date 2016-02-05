function Gamma = update_closest_simp(Gamma,Surface)

J=size(Gamma.X,1);
Nsimp = size(Surface.tri,1); 

for i=1:J
    
    
    % Consider former closest simplex
    i_simp = Gamma.closest_simp(i,1); 
    
    [inside,q,delta_lambda] = is_inside(Gamma.X(i,:),i_simp,Surface);
    
    if(inside)
        % Use projection node
        Gamma.X(i,:) = q;

    else
        % Save delta_lambda and projected node 
        delta_lambda_min = delta_lambda;
        q_min = q;
        i_min = i_simp;
        
        % Look for neighbors
        used = zeros(Nsimp,1);
        found = 0;
        found_q = zeros(3,1); 
        
        level = 1;
    
        list = Surface.neigh(i_simp,:)';
        
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
                [inside,q,delta_lambda] = is_inside(Gamma.X(i,:),i_simp_list,Surface);
                    
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
            Gamma.closest_simp(i,1) = found;
            Gamma.X(i,:) = found_q;
            

        else
            fprintf('\nCannot find a closest simplex for X(%d) = (%3.5f, %3.5f, %3.5f)\n', i, Gamma.X(i,1),Gamma.X(i,2),Gamma.X(i,3)); 
            pause
        end
    end
    
end

% Update for Lambda
for k=1:size(Gamma.Lambda,1)
    p = Gamma.X(Gamma.Lambda(k,1),:); 
    isimp = Gamma.closest_simp(Gamma.Lambda(k,1),1); 
    
    for l=2:3
        Gamma.X(Gamma.Lambda(k,l),:) = p;
        Gamma.closest_simp(Gamma.Lambda(k,l),1) = isimp;
    end
end


end




