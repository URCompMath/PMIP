function [Gamma_new,S,P,M,b] = solve_umfpack(Gamma,b,delta_t,config, Image)


[P,M]  = get_sparse_projection_matrix(Gamma,Image);
S  = get_sparse_system_matrix(Gamma,config,delta_t);


% S = P*S*P; 
% b = P*b;
% X = S\b;


% Solution of PSPx = Pb, Solution x lies in subspace {z = Pz}
% (1) Solve Sy = b, solution y fullfills automatically PSY = Pb
% (2) Set x = Py --> Px = PPy = Py = x, as (PP=P)

S = P*S*P; 
b = P*b; 
X = (S+M)\b;


Gamma_new = Gamma; 
Gamma_new.delta_X = reshape(X,2,size(Gamma.X_old,1))';

end
