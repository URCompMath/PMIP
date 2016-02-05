R = 0.9;
r = 0.5;
for i=1:size(Gamma.X,1)
    x = Gamma.X(i,:);
    v = [x(1),x(2),0];
    x0 = R*v/norm(v);
    y = x - x0;
    xnew = x0 + r*y/norm(y);
    Gamma.X(i,:) = xnew;
end
