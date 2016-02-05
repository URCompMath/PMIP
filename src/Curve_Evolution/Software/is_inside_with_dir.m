function [in,q,delta_lambda] = is_inside_with_dir(p,v,i_simp,Surface)

p1 = Surface.nodes(Surface.tri(i_simp,1),:);
p2 = Surface.nodes(Surface.tri(i_simp,2),:);
p3 = Surface.nodes(Surface.tri(i_simp,3),:);


% Intersection:  line p + lambda * v with simplex plane
v1 = p2-p1;
v1 = v1/norm(v1);
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
delta_lambda = 0;

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
        
        if(lambda(j)>1)
            delta_lambda = delta_lambda + (lambda(j)-1)^2;
        else
            delta_lambda = delta_lambda + (lambda(j))^2;
        end
        
    end
end


end