function boundary_index=is_boundary(x,nx,ny)

tol = 1.5; 

p = [x(1),0.5; nx+0.5,x(2); x(1),ny+0.5; 0.5,x(2)];


diff = ones(4,1)*x-p;


[val,k]=min(sum(diff.*diff,2)); 
if(val < tol^2)
    boundary_index = k; 
else
    boundary_index = 0;
end
end