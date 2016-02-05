function [S,b] = get_sparse_system_matrix_and_right_hand_side(config,Gamma,Omega,compute_external)

% Number of nodes, number of elements
J = size(Gamma.X,1); 
J0 = size(Gamma.simplices,1); 
 
%% Step 1: Preparation: Average of neighbor centers
% For config.mesh.flag == 2 (BGN strategy) compute average of neighbors
if(config.mesh.flag == 2)
    z_avg = zeros(J,3); 
    n_z   = zeros(J,1); 
    
    if(config.mesh.lambda_t>0)
        % Average of simplex centers
        for i=1:size(Gamma.simplices,1)
            i1 = Gamma.simplices{i,1}.nodes(1);
            i2 = Gamma.simplices{i,1}.nodes(2);
            i3 = Gamma.simplices{i,1}.nodes(3);
            
            c = (Gamma.X(i1,:) + Gamma.X(i2,:) + Gamma.X(i3,:))/3;
            z_avg(i1,:) = z_avg(i1,:) + c;
            z_avg(i2,:) = z_avg(i2,:) + c;
            z_avg(i3,:) = z_avg(i3,:) + c;
            
            n_z(i1) = n_z(i1)+1;
            n_z(i2) = n_z(i2)+1;
            n_z(i3) = n_z(i3)+1;
            
            
        end
        
        % Divide by n_z
        for i=1:J
            z_avg(i,:) = z_avg(i,:) / n_z(i);
        end
    end
end    
 
%% Step 2: Calculate b, c1, c2
b  = zeros(J,1); 
c1 = zeros(J,1);
c2 = zeros(J,1); 
 
if(compute_external)
    for j=1:J0
        nodes = Gamma.simplices{j,1}.nodes;
        for k=1:3
            b(nodes(k)) = b(nodes(k)) + Gamma.area_sigma(j)/3*func_f(Gamma.X(nodes(k),:)',Gamma,Omega,config, Gamma.simplices{j,1}.index(1));
            
            if(config.mesh.flag == 2 && max(config.mesh.sigma_t)>0 && config.mesh.lambda_t > 0)  % BGN strategy
                c1(nodes(k)) = c1(nodes(k)) + Gamma.area_sigma(j)/3 * config.mesh.lambda_t /Gamma.delta_t * (z_avg(nodes(k),:) - Gamma.X(nodes(k),:)) * (Gamma.tau1(nodes(k),:)');
                c2(nodes(k)) = c2(nodes(k)) + Gamma.area_sigma(j)/3 * config.mesh.lambda_t /Gamma.delta_t * (z_avg(nodes(k),:) - Gamma.X(nodes(k),:)) * (Gamma.tau2(nodes(k),:)');
            end
        end
    end
end

 
%% Step 3: Calculate 
% Calculate 1/sigma M^{-1}b
help = 1/config.method.sigma * (3./ Gamma.area_Lambda).*b; 
 
 
% Calculate 1/sigma_t M^{-1}ci
helpt1 = zeros(J,1);
helpt2 = zeros(J,1);
if(config.mesh.flag ==2 && max(config.mesh.sigma_t)>0)
    for j=1:J
        if(config.mesh.sigma_t(j,1) > 0)
            helpt1(j,1) =  ((3/Gamma.area_Lambda(j,1))*c1(j,1) ) / config.mesh.sigma_t(j,1);
            helpt2(j,1) =  ((3/Gamma.area_Lambda(j,1))*c2(j,1) ) / config.mesh.sigma_t(j,1);
        end
    end
end

% Prepare AX_old and N*help, Tihelpti (i=1,2)
% A*X_old
Ax = zeros(J,3); 
 
NMb = zeros(J,3); 
T1Mc1 = zeros(J,3); 
T2Mc2 = zeros(J,3); 
 
 
% Prepare NM^{-1}N^T, T_i M^{-1} T_i^T (i=1,2)
row = zeros(9*J,1); 
col = zeros(9*J,1); 
 
for j=1:J
    row((9*j-8):(9*j),1) = [3*j-2; 3*j-2; 3*j-2; 3*j-1; 3*j-1; 3*j-1; 3*j; 3*j; 3*j]; 
    col((9*j-8):(9*j),1) = [3*j-2; 3*j-1; 3*j; 3*j-2; 3*j-1; 3*j; 3*j-2; 3*j-1; 3*j];
end

 
val_NMN = zeros(9*J,1); 
val_T1MT1 = zeros(9*J,1); 
val_T2MT2 = zeros(9*J,1); 
 
% Calculate A
row_A = zeros(24*3*J,1); 
col_A = zeros(24*3*J,1); 
val_A = zeros(24*3*J,1); 
count_A = 0; 
 
 
 
% Loop over all simplices, compute contribution of each simplex to right
% hand side and system matrix
for j=1:J0

    nodes = Gamma.simplices{j,1}.nodes;
    nabla_sigma_j = reshape(Gamma.nabla_sigma(j,:),3,3)'; 
    
    % Add contribution of sigma_j to vectors and matrices
    for k=1:3
        for l=1:3       
            % Compute AXold
            Ax(nodes(k),:)= Ax(nodes(k),:) + Gamma.area_sigma(j)* (nabla_sigma_j(k,:)*nabla_sigma_j(l,:)')*Gamma.X(nodes(l),:); 
            
            % Compute A
            row_loc = [3*nodes(k)-2; 3*nodes(k)-1; 3*nodes(k)]; 
            col_loc = [3*nodes(l)-2; 3*nodes(l)-1; 3*nodes(l)];
            val_loc = Gamma.area_sigma(j)* nabla_sigma_j(k,:)*(nabla_sigma_j(l,:)') * ones(3,1); 
            
            count_A = count_A + 3; 
            row_A((count_A-2):count_A,1) = row_loc; 
            col_A((count_A-2):count_A,1) = col_loc;
            val_A((count_A-2):count_A,1) = val_loc; 
        end

        % Compute NMb, T1Mc1, T2Mc2 
        NMb(nodes(k),:)   = NMb(nodes(k),:)   + 1/3*Gamma.area_sigma(j)*Gamma.nu_sigma(j,:)*help(nodes(k)); 
        
        if(config.mesh.flag == 2 && max(config.mesh.sigma_t)>0)
            T1Mc1(nodes(k),:) = T1Mc1(nodes(k),:) + 1/3*Gamma.area_sigma(j)*Gamma.tau1(nodes(k),:) * helpt1(nodes(k));
            T2Mc2(nodes(k),:) = T2Mc2(nodes(k),:) + 1/3*Gamma.area_sigma(j)*Gamma.tau2(nodes(k),:) * helpt2(nodes(k));
        end
        
        % Compute NM^{-1}N^T, TiM^{-1}Ti^T, i=1,2
        add_mat = 1/3*Gamma.area_sigma(j) * (Gamma.omega(nodes(k),:)')*(Gamma.omega(nodes(k),:)); 
        add_vec = reshape(add_mat',9,1); 
        val_NMN((9*nodes(k)-8):(9*nodes(k)),1) = val_NMN((9*nodes(k)-8):(9*nodes(k)),1) + add_vec; 
        
        if(config.mesh.flag == 2 && max(config.mesh.sigma_t)>0)
            add_mat1 = 1/3*Gamma.area_sigma(j) * (Gamma.tau1(nodes(k),:)')*(Gamma.tau1(nodes(k),:));
            add_mat2 = 1/3*Gamma.area_sigma(j) * (Gamma.tau2(nodes(k),:)')*(Gamma.tau2(nodes(k),:));
            
            if(config.mesh.sigma_t(nodes(k),1) > 0)
                add_vec1 = reshape(add_mat1',9,1) / config.mesh.sigma_t(nodes(k),1); 
                add_vec2 = reshape(add_mat2',9,1) / config.mesh.sigma_t(nodes(k),1); 
            else
                add_vec1 = zeros(9,1); 
                add_vec2 = zeros(9,1); 
            end
            val_T1MT1((9*nodes(k)-8):(9*nodes(k)),1) = val_T1MT1((9*nodes(k)-8):(9*nodes(k)),1) + add_vec1;
            val_T2MT2((9*nodes(k)-8):(9*nodes(k)),1) = val_T2MT2((9*nodes(k)-8):(9*nodes(k)),1) + add_vec2;
            
            
        end
        
        
        
    end
end

 
% Convert right hand side (output) vector
b = zeros(3*J,1); 
 
for i=1:J
    b(3*i-2,1)=-Ax(i,1) + NMb(i,1) + T1Mc1(i,1) + T2Mc2(i,1);  % if config.mesh.flag \neq 2, T1Mc1 and T2Mc2 store 0! 
    b(3*i-1,1)=-Ax(i,2) + NMb(i,2) + T1Mc1(i,2) + T2Mc2(i,2);
    b(3*i,1)  =-Ax(i,3) + NMb(i,3) + T1Mc1(i,3) + T2Mc2(i,3);
end

% Convert system matrix (output)
row_A = row_A(1:count_A,1); 
col_A = col_A(1:count_A,1); 
val_A = val_A(1:count_A,1); 
 
A     = sparse(row_A,col_A,val_A, 3*J, 3*J); 
NMN   = sparse(row,col,val_NMN, 3*J, 3*J); 
T1MT1 = sparse(row,col,val_T1MT1, 3*J, 3*J); 
T2MT2 = sparse(row,col,val_T2MT2, 3*J, 3*J); 
 
S = 1/Gamma.delta_t * (1/config.method.sigma * NMN + T1MT1 + T2MT2) + A;  
 
 
end
