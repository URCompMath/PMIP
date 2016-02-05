R = 1;
for i=1:size(Gamma.X,1)
    x = Gamma.X(i,:);
    Gamma.X(i,:) = R*x/norm(x);
end