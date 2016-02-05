function [Gamma_new,N_nodes_new] = refine_curve(Gamma,Image,Surface,node_i0,N_nodes)
% initialize output
Gamma_new = Gamma; 
N_nodes_new = N_nodes; 

% Start loop over nodes belonging to sub-curve with representive node_i0
% Add a new point in the middle of two neighboring points, project to
% surface in case dimension==3
i = node_i0; 
j = 1; 

J = size(Gamma.X,1); 

while(i~=node_i0 || j==1)
    if(Gamma.neigh(i,2)>0)
        xi  = Gamma.X(i,:);
        xii = Gamma.X(Gamma.neigh(i,2),:); 
        % calc new point in the middle
        m = (xi+xii)/2;
        % project to surface if dimension_flag = 3
        if(size(Gamma.X,2) == 3)
            [i_simp_closest,m] = get_closest_simp_use_parents(m,Gamma.closest_simp(i,1),Surface);
        end
        J = J + 1; 
        % Add new point to curve 
        Gamma_new.X(J,:) = m;
        % Adapt neighbor information
        Gamma_new.neigh(i,2) = J;
        Gamma_new.neigh(Gamma.neigh(i,2),1) = J;
        % Add neigh info
        Gamma_new.neigh(J,1) = i;
        Gamma_new.neigh(J,2) = Gamma.neigh(i,2); 
        
        if(size(Gamma.X,2)==3)
            % Add closest_simp entry
            Gamma_new.closest_simp(J,1) = i_simp_closest;
        end
        
        N_nodes_new = N_nodes_new + 1;

        % Go to next point
        i = Gamma.neigh(i,2);
    else
        i = node_i0; % end of curve --> stop loop
    end
    
    j = j+1;
end



end