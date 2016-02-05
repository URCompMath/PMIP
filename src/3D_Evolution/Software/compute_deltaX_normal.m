function dXn = compute_deltaX_normal(Gamma)
J = size(Gamma.delta_X,1); 

dXn = zeros(J,1); 
for i=1:J
    dXn(i,1) = Gamma.delta_X(i,:) * Gamma.omega(i,:)'; 
end
end