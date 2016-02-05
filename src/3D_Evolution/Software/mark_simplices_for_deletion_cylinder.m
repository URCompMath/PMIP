function mark = mark_simplices_for_deletion_cylinder(X,simplices,p0,n0,a)

N = size(simplices,1); 
h = 10*a;
r = 10*a;

mark = zeros(N,1); 

n0=n0/norm(n0);

for i=1:N
    inside = 0; 
    for j=1:3
        x = X(simplices{i,1}.nodes(j),:);
        
        v = x - p0'; % note: p0: 3x1 vector, n0: 3x1 vec
        vn = v*n0; 
        if(abs(vn) <= h)
            vt = v - vn*(n0'); 
            if(norm(vt) <= r)
                inside = 1;
            end
        end
    end
    if(inside == 1)
        mark(i) = 1;
    end
    
end

end