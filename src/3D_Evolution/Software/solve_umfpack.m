function Gamma_new = solve_umfpack(config,Gamma,Omega,compute_external)

% Get system matrix (sparse) and right hand side
[S,b] = get_sparse_system_matrix_and_right_hand_side(config,Gamma,Omega,compute_external); 

% Solve Linear Equation using Matlab built-in routine
X = S\b;

% Output generation
Gamma_new = Gamma; 
Gamma_new.delta_X = reshape(X,3,size(Gamma.X,1))'; 
end
