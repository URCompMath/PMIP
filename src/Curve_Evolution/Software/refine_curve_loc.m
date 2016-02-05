function Gamma_new = refine_curve_loc(Gamma,Surface,i0)
% Insert a point between i0 and neigh(i0)

% initialize output
Gamma_new = Gamma;

% Get position corresponding to i0 and neigh(i0,2)
xi = Gamma.X(i0,:);
xii = Gamma.X(Gamma.neigh(i0,2),:);

J = size(Gamma.X,1);

% calc new point in the middle
m = (xi+xii)/2;
% project to surface if dimension_flag = 3
if(size(Gamma.X,2) == 3)
      [i_simp_closest,m] = get_closest_simp_use_parents(m,Gamma.closest_simp(i0,1),Surface); 
end
J = J + 1;

% Add new point to curve
Gamma_new.X(J,:) = m;
% Adapt neighbor information
Gamma_new.neigh(i0,2) = J;
Gamma_new.neigh(Gamma.neigh(i0,2),1) = J;
% Add neigh info
Gamma_new.neigh(J,1) = i0;
Gamma_new.neigh(J,2) = Gamma.neigh(i0,2);

if(size(Gamma.X,2)==3)
    % Add closest_simp entry
    Gamma_new.closest_simp(J,1) = i_simp_closest;
end

end