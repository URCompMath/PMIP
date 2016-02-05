function [n,p0] = get_separation_plane(X,index)

p0 = zeros(3,1); 
std0 = zeros(3,1); 
for i=1:3
    std0(i)=std(X(index,i));
    p0(i)  =mean(X(index,i)); 
end

[std_min,k] = min(std0); 
switch k
    case 1
        perm_123 = [2;3;1]; 
    case 2
        perm_123 = [3;1;2];
    case 3
        perm_123 = [1;2;3];
end

N = size(index,1); 
A = zeros(3,3); 
y = zeros(3,1); 

for i=1:N
    A(1,1) = A(1,1) + X(index(i),perm_123(1))^2;
    A(1,2) = A(1,2) + X(index(i),perm_123(1))*X(index(i),perm_123(2));
    A(1,3) = A(1,3) + X(index(i),perm_123(1)); 
    A(2,2) = A(2,2) + X(index(i),perm_123(2))^2; 
    A(2,3) = A(2,3) + X(index(i),perm_123(2));
    y(1)   = y(1)   + X(index(i),perm_123(1))*X(index(i),perm_123(3));
    y(2)   = y(2)   + X(index(i),perm_123(2))*X(index(i),perm_123(3));
    y(3)   = y(3)   + X(index(i),perm_123(3));
    
end
A(2,1)=A(1,2);
A(3,1)=A(1,3);
A(3,2)=A(2,3);
A(3,3)=N; 

a = A\y;

p0(perm_123(3)) = a(1)*p0(perm_123(1)) + a(2)*p0(perm_123(2)) + a(3); 

switch k
    case 1
        n = [1;-a(1);-a(2)];
    case 2
        n = [-a(2);1;-a(1)];
    case 3
        n = [-a(1);-a(2);1];
end
n = n/norm(n); 
end



