function u = diffusion_eq_geodesic2(Surface,Image,Omega,lambda)
 
% 1/lambda < nabla_M u, nabla_M phi > + < u, phi > = < u0, phi > for all
% phi in S^h

% Preparation 
[area_sigma, nabla_sigma] = compute_help_quantities(Surface); 
Surface.area_sigma = area_sigma;
Surface.nabla_sigma = nabla_sigma; 

NR = size(Omega.coeffs,1); 
Nsimp = size(Surface.tri,1); % Image approximation
u = zeros(size(Image.data));

color_flag = 0; 
if(size(Image.data,2)==3)
    color_flag = 1;
end

for k=1:NR
    fprintf('Region %d/%d\n', k,NR); 
    % Compute right hand side
    [b,inside] = compute_right_hand_side_diffusion_geodesic(Surface,Image,Omega,k);
    
    % Compute sparse system matrix
    S = get_sparse_system_matrix_diffusion_geodesic(Surface,Omega,lambda,inside,k);
    
    % Solve diffusion equation with finite elements, compute image data on the
    % nodes
    if(color_flag == 0)
        U = S\b;
    else
        U1 = S\b(:,1);
        U2 = S\b(:,2);
        U3 = S\b(:,3); 
        
        U = [U1,U2,U3];
    end
    
    % Get image approximation on the simplices
    for i=1:Nsimp
        if(Omega.A_info(i,1)==k)
            nodes = Surface.tri(i,:);
            u(i,:) = (U(nodes(1),:) + U(nodes(2),:) + U(nodes(3),:))/3;
        end
    end
end

 
end

function [area_sigma, nabla_sigma] = compute_help_quantities(Surface)
Nsimp = size(Surface.tri,1); 
 
area_sigma = zeros(Nsimp,1); 
nabla_sigma = zeros(Nsimp,9); 
 
for j=1:Nsimp
    % Compute area of element sigma_j
    nodes = Surface.tri(j,:); 
    X1 = Surface.nodes(nodes(1),:); 
    X2 = Surface.nodes(nodes(2),:); 
    X3 = Surface.nodes(nodes(3),:); 
    
    cross_j = cross(X1-X3,X2-X3); 
    area_sigma(j,1) = 0.5*norm(cross_j); 
    
    % Compute nabla_sigma
    % for chi_1 let tau_1 = (X3-X2)/|X3-X2|, 
    %               tau_2 = (X1-q)/|X1-q|, q = X2 + (X1-X2)*tau_1 tau_1
    
    % nabla_tau_1 chi_1 = 0
    % nabla_tau_2 chi_1 = 1/(X1-q); 
    % nabla_M chi_1 = nabla_tau_2 chi_1 * tau_2 = (X1-q)/|X1-q|^2
    
    t1 = X3-X2;
    t1 = t1/norm(t1); 
    q = X2 + ((X1-X2)*(t1'))*t1;
    u = X1-q;
    nabla_sigma(j,1:3) = u/(norm(u)^2); 
    
    t1 = X1-X3;
    t1 = t1/norm(t1); 
    q = X3 + ((X2-X3)*(t1'))*t1;
    u = X2-q;
    nabla_sigma(j,4:6) = u/(norm(u)^2); 
    
    t1 = X2-X1;
    t1 = t1/norm(t1); 
    q = X1 + ((X3-X1)*(t1'))*t1;
    u = X3-q;
    nabla_sigma(j,7:9) = u/(norm(u)^2); 
    

end

end

 
function [b,inside] = compute_right_hand_side_diffusion_geodesic(Surface,Image,Omega,region_index)
% b_j = < chi_j, u_0 >

 
Nnodes = size(Surface.nodes,1); 
Nsimp  = size(Surface.tri,1); 

dim = size(Image.data,2); 

b = zeros(Nnodes,dim); 
inside = zeros(Nnodes,1); 
 
for i=1:Nsimp
    % Compute contribution of sigma_i 
    if(Omega.A_info(i,1)==region_index)
        nodes = Surface.tri(i,:);
        for j=1:3
            b(nodes(j),:) = b(nodes(j),:) + Surface.area_sigma(i,1) * 1/3 * Image.data(i,:);
            inside(nodes(j),1) = 1; 
        end
    end
end

end


function S = get_sparse_system_matrix_diffusion_geodesic(Surface,Omega,lambda,inside,region_index)

Nnodes = size(Surface.nodes,1); 
Nsimp  = size(Surface.tri,1); 

row = (1:Nnodes)'; 
col = row;
val = ones(Nnodes,1)-inside;  % 1 at nodes outside the region --> eq: x_i = 0 since b_i = 0 outside the region 

N = Nnodes;


for i=1:Nsimp
    if(Omega.A_info(i,1) == region_index)
        nodes = Surface.tri(i,:); 
        
        nabla_sigma_i = reshape(Surface.nabla_sigma(i,:), 3, 3)';
        for j=1:3
            val(nodes(j),1) = val(nodes(j),1) + 1/(3) * Surface.area_sigma(i,1);
            for l=1:3
                N=N+1;
                row(N,1) = nodes(j); 
                col(N,1) = nodes(l); 
                val(N,1) = Surface.area_sigma(i,1) * 1/lambda *(nabla_sigma_i(j,:)*nabla_sigma_i(l,:)');
            end
        end
    end
end
S = sparse(row,col,val,Nnodes,Nnodes); 
        


end
 



