function [cp_out,i_simp] = project2surface_init(p,v,Surface)

N = 0;
list = zeros(20,4); % (i_simp, q(1),q(2),q(3))

Nsimp = size(Surface.tri,1); 



for i=1:Nsimp
    % Intersection:  line p + lambda * v with simplex plane
    p1 = Surface.nodes(Surface.tri(i,1),:) ;
    p2 = Surface.nodes(Surface.tri(i,2),:) ;
    p3 = Surface.nodes(Surface.tri(i,3),:) ;
    
    v1 = p2-p1;
    v1 = v1/norm(v1) ;
    v2 = p3-p1; 
    v2 = v2/norm(v2);
       
    
    M = [-v',v1',v2'];
    l = M \ ((p-p1)');
    
    q = p + l(1)*v;
    
    % Compute barycentric coordinates of q
    b = [q(1)-p3(1); q(2)-p3(2)];
    A = [p1(1)-p3(1), p2(1)-p3(1); p1(2)-p3(2), p2(2)-p3(2)];
    if(abs(det(A))<1e-10)
        b = [q(1)-p3(1); q(3)-p3(3)];
        A = [p1(1)-p3(1), p2(1)-p3(1); p1(3)-p3(3), p2(3)-p3(3)];
        if(abs(det(A))<1e-10)
            b = [q(2)-p3(2); q(3)-p3(3)];
            A = [p1(2)-p3(2), p2(2)-p3(2); p1(3)-p3(3), p2(3)-p3(3)];
        end
    end
    
    
    lambda = A\b;
    lambda(3)=1-lambda(2)-lambda(1);
    
    
    in = 1;
    for j=1:3
        % round using tolerance
        if(lambda(j)>1 && lambda(j)<=1+Surface.tol)
            lambda(j)=1;
        else
            if(lambda(j)<0 && lambda(j)>= -Surface.tol)
                lambda(j)=0;
            end
        end
    
        % Set in to 0 if point is outside
        if(lambda(j)>1 || lambda(j)<0)
            in = 0;
        end
    end
    
    
    if(in==1)
        % Projection q is inside the simplex, possible candidate
        N=N+1;
        list(N,:) = [i, q]; 
    end
end

% Find candidate q from list(:,2:4) which is the closest
dist_vec = list(1:N,2:4) - ones(N,1) * p;
dist = sqrt(sum(dist_vec .* dist_vec,2)); 
[val,ind] = min(dist); 

cp_out = list(ind,2:4);
i_simp = list(ind,1);
end