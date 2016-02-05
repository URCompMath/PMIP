function [cp,i_simp] = get_closest_point_and_simplex(p,Surface)
% Get closest point on surface and simplex of triangulated surface
% Function used by func_I



ap = zeros(1,3);
ap(1) = (p(1)-Surface.aminmax(1))/Surface.grid_width + 1;
ap(2) = (p(2)-Surface.aminmax(3))/Surface.grid_width + 1;
ap(3) = (p(3)-Surface.aminmax(5))/Surface.grid_width + 1;

a0 = zeros(1,3);
a1 = zeros(1,3);
for j=1:3
    a0(j)=floor(ap(j));
    a1(j)=ceil(ap(j));
end

nsimp = 0;
width = 0;

while(nsimp==0 && width <= 2)
    Nx0 = max([1,a0(1)-width]);
    Nx1 = min([a1(1)+width, size(Surface.grid_info,1)]);
    Ny0 = max([1,a0(2)-width]);
    Ny1 = min([a1(2)+width, size(Surface.grid_info,2)]);
    Nz0 = max([1,a0(3)-width]);
    Nz1 = min([a1(3)+width, size(Surface.grid_info,3)]);
    
    Nc = (Nx1-Nx0+1)*(Ny1-Ny0+1)*(Nz1-Nz0+1);
    corners = zeros(Nc,3);
    count = 1;
    
    for i=Nx0:Nx1
        for j=Ny0:Ny1
            for k=Nz0:Nz1
                corners(count,:) = [i,j,k];
                count = count + 1;
            end
        end
    end
    
    % corners
    
    
    nsimp = 0;
    list = zeros(0,1);
    
    % Collect simplices
    for j=1:Nc
        q = corners(j,:);
        Surface.grid_info{q(1), q(2), q(3)};
        if(~isempty(Surface.grid_info{q(1), q(2), q(3)}))
            nsimp_j = Surface.grid_info{q(1), q(2), q(3)}.nsimp;
            index_simplex_j = Surface.grid_info{q(1), q(2), q(3)}.index_simplex(1:nsimp_j,1);
            
            for k=1:nsimp_j
                index_jk = index_simplex_j(k,1);
                test_i = find(list==index_jk);
                if(size(test_i,1)==0)
                    nsimp = nsimp+1;
                    list(nsimp,1)=index_jk;
                end
            end
        end
    end
    width = width + 1;
end


% Find closest simplex i from list
i_simp = -1;
dist = 1e8;
cp = zeros(1,3);

for i=1:nsimp
    p1 = Surface.nodes(Surface.tri(list(i,1),1),:);
    p2 = Surface.nodes(Surface.tri(list(i,1),2),:);
    p3 = Surface.nodes(Surface.tri(list(i,1),3),:);
    
    [dist_i,cp_i] = dist2simplex(p,p1,p2,p3);
    if(dist_i<dist)
        dist = dist_i;
        i_simp = list(i,1);
        cp = cp_i;
    end
end
end


