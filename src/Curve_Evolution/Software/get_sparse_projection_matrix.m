function [P,M] = get_sparse_projection_matrix(Gamma,Image)
J = size(Gamma.X,1);

% Sparse Matrix P
row = zeros(4*J,1);
col = zeros(4*J,1);
val = zeros(4*J,1);


rowM = zeros(2*J,1); 
colM = zeros(2*J,1); 
valM = zeros(2*J,1); 
countM = 0; 

% Prepare mat_n
n1 = [0;-1];
n2 = [1;0];
n3 = [0;1];
n4 = [-1;0];
n = [n1,n2,n3,n4];

mat_n = zeros(4,4);
for i=1:4
    help = n(:,i)*n(:,i)';
    mat_n(:,i) = [help(1,1);help(1,2);help(2,1);help(2,2)];
end



% Set P
for i=1:J
    row(((4*i-3):(4*i))') = [2*i-1; 2*i-1; 2*i; 2*i];
    col(((4*i-3):(4*i))') = [2*i-1; 2*i; 2*i-1; 2*i];
    
    % Boundary points
    if(Gamma.neigh(i,1) == -1 || Gamma.neigh(i,2)==-1)
        ind = is_boundary(Gamma.X_old(i,:),Image.sizes(1),Image.sizes(2));
        val(((4*i-3):(4*i))') = [1;0;0;1] - mat_n(:,ind);
        
        n_omega = n(:,ind); 
        
        if(abs(n_omega(1))>1e-4)
            rowM(((countM+1):(countM+2))') = [2*i-1; 2*i-1];
            colM(((countM+1):(countM+2))') = [2*i-1; 2*i];
            valM(((countM+1):(countM+2))') = [1; n_omega(2)/n_omega(1)];
        else
            rowM(((countM+1):(countM+2))') = [2*i; 2*i];
            colM(((countM+1):(countM+2))') = [2*i-1; 2*i];
            valM(((countM+1):(countM+2))') = [n_omega(1)/n_omega(2);1];
        end
        countM=countM+2; 
           
        
    else
        val(((4*i-3):(4*i))') = [1;0;0;1];
    end
end
rowM = rowM(1:countM,1); 
colM = colM(1:countM,1); 
valM = valM(1:countM,1); 





P = sparse(row,col,val,2*J,2*J);
M = sparse(rowM,colM,valM,2*J,2*J); 

end