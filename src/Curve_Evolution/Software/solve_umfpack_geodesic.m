function Gamma_new = solve_umfpack_geodesic(Gamma,b,delta_t,config,repeat)

% Step 1: Compute sparse system matrix S
% S = get_sparse_system_matrix_geodesic(Gamma,config,delta_t);
[S,P] = get_sparse_system_matrix_geodesic2(Gamma,config,delta_t);
    
% Apply P on b
b = P*b;


% % Step 2: Solve linear system with UMFPACK (Matlab built-in)
%X = S \ b; 

% Alternativ: Solve with cg since S is singular (sub-space)
X=solve_geodesic_cg(S,b,config,repeat);


J = size(Gamma.X,1); 
Gamma.delta_X = (reshape(X,3,J))';


Gamma_new = Gamma;
end


