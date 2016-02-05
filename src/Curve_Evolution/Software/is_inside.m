function [in,q,delta_lambda] = is_inside(p,i_simp,Surface)

p1 = Surface.nodes(Surface.tri(i_simp,1),:);
p2 = Surface.nodes(Surface.tri(i_simp,2),:);
p3 = Surface.nodes(Surface.tri(i_simp,3),:);


% Project p to simplex plane
n1 = p2-p1;
n1 = n1/norm(n1);
n2 = p3-p1;
n2 = n2/norm(n2);

n = cross(n1,n2);

q = p - ((p-p1)*n')*n;

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