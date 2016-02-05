function nu_F = calc_normal_to_surface(i_node,Gamma,Surface)
% Assumes that all vertices of an triangle are ordered anti-clockwise 
% Compute nu_F, the outer unit normal to the surface


i_simp = Gamma.closest_simp(i_node,1); 

p1 = Surface.nodes(Surface.tri(i_simp,1),:);
p2 = Surface.nodes(Surface.tri(i_simp,2),:);
p3 = Surface.nodes(Surface.tri(i_simp,3),:); 

u = p2-p1; 
u = u/norm(u); 

v = p3-p1; 
v = v/norm(v); 

nu_F = cross(u,v); 
nu_F = nu_F/norm(nu_F); 
end