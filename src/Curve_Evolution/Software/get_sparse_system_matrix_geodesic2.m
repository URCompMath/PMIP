function [S,P] = get_sparse_system_matrix_geodesic2(Gamma,config,delta_t)

J=size(Gamma.X,1); 

%% Preparation for triple junctions
is_triple = zeros(J,3); 
N_triple = size(Gamma.Lambda,1); 
ind_nr = [2,3; 3,1; 1,2]; 

for i=1:N_triple
    nodes = Gamma.Lambda(i,1:3); 
    for j=1:3
        is_triple(nodes(j),:) = [1, nodes(ind_nr(j,1)), nodes(ind_nr(j,2))];
    end
end
    

%% Get Sparse Matrices 

row_P_geo = zeros(9*J,1); 
col_P_geo = zeros(9*J,1); 
val_P_geo = zeros(9*J,1); 

row_P_triple = zeros(3*J,1); 
col_P_triple = zeros(3*J,1); 
val_P_triple = zeros(3*J,1); 

row = zeros(3*J,1); 
col = zeros(3*J,1);
val = zeros(3*J,1); 

ig = 0;
it = 0; 
i0 = 0; 

for i=1:J
    % 1) P_geo = (R^3)^J --> { (X1,...,XJ) : Xi . nu_F = 0}    
    row_P_geo(((ig+1):(ig+9))',1) = [3*i-2; 3*i-2; 3*i-2; 3*i-1; 3*i-1; 3*i-1; 3*i  ; 3*i  ; 3*i  ];
    col_P_geo(((ig+1):(ig+9))',1) = [3*i-2; 3*i-1; 3*i;   3*i-2; 3*i-1; 3*i;   3*i-2; 3*i-1; 3*i  ];
       
    nuF = Gamma.nu_F(i,:); 
    M_block = eye(3,3) - nuF'*nuF;
    val_P_geo(((ig+1):(ig+9))',1) = reshape(M_block',9,1); 

    ig = ig + 9;

    % 2) P_triple = (R^3)^J --> {(X1,...,XJ) : X_ik1 = X_ik2 = X_ik3 for all k=1,...,NT}
    row_P_triple(((it+1):(it+3))',1) = [3*i-2; 3*i-1; 3*i  ];
    col_P_triple(((it+1):(it+3))',1) = [3*i-2; 3*i-1; 3*i  ];

    if(is_triple(i,1) == 0)
        % nodes does not belong to a triple junction
        val_P_triple(((it+1):(it+3))',1) = [1; 1; 1]; 
    else
        % nodes belongs to a triple junction
        val_P_triple(((it+1):(it+3))',1) = 1/3 * [1; 1; 1]; 
        
        % Two further 1/3 eye(3,3) blocks
        i2 = is_triple(i,2); 
        i3 = is_triple(i,3); 
        
        row_P_triple(((it+4):(it+6))',1) = [3*i-2;  3*i-1;  3*i  ];
        col_P_triple(((it+4):(it+6))',1) = [3*i2-2; 3*i2-1; 3*i2  ];
        val_P_triple(((it+4):(it+6))',1) = 1/3 * [1; 1; 1]; 
        
        row_P_triple(((it+7):(it+9))',1) = [3*i-2;  3*i-1;  3*i  ];
        col_P_triple(((it+7):(it+9))',1) = [3*i3-2; 3*i3-1; 3*i3  ];
        val_P_triple(((it+7):(it+9))',1) = 1/3 * [1; 1; 1]; 
        
        it = it + 6;
    end
    
    it = it + 3;
    
    % 3) S = 1/sigma tau_m N_M M^{-1} N_M + A
    omega = Gamma.nu(i,:); 
    
    if(Gamma.neigh(i,1) > 0 && Gamma.neigh(i,2) > 0)
        % inner point
        h1 = norm(Gamma.X_old(i,:)-Gamma.X_old(Gamma.neigh(i,1),:));
        h2 = norm(Gamma.X_old(Gamma.neigh(i,2),:)-Gamma.X_old(i,:)); 
        
        factor = (h1+h2)/2;
        
        % S_ii
        row(((i0+1):(i0+9))') = [3*i-2; 3*i-2; 3*i-2; 3*i-1; 3*i-1; 3*i-1; 3*i  ; 3*i  ; 3*i  ];
        col(((i0+1):(i0+9))') = [3*i-2; 3*i-1; 3*i  ; 3*i-2; 3*i-1; 3*i  ; 3*i-2; 3*i-1; 3*i  ];
        
        M_block = 1/(config.method.sigma*delta_t)*factor* ((omega')*omega) + (1/h1 + 1/h2) * eye(3,3);
        val(((i0+1):(i0+9))') = reshape(M_block',9,1);
        
        % S_iim and S_iip
        im = Gamma.neigh(i,1); 
        ip = Gamma.neigh(i,2); 
        
        row(((i0+10):(i0+12))') = [3*i-2;  3*i-1;  3*i ];
        col(((i0+10):(i0+12))') = [3*im-2; 3*im-1; 3*im];
        val(((i0+10):(i0+12))') = -1/h1 * [1; 1; 1]; 
        
        row(((i0+13):(i0+15))') = [3*i-2;  3*i-1;  3*i ];
        col(((i0+13):(i0+15))') = [3*ip-2; 3*ip-1; 3*ip];
        val(((i0+13):(i0+15))') = -1/h2 * [1; 1; 1]; 
        
        i0 = i0 + 15; 
        
    else
        if(Gamma.neigh(i,1) < 0)
            h = norm(Gamma.X_old(Gamma.neigh(i,2),:)-Gamma.X_old(i,:)); 
            factor = h/2; 
            
            % S_ii 
            row(((i0+1):(i0+9))') = [3*i-2; 3*i-2; 3*i-2; 3*i-1; 3*i-1; 3*i-1; 3*i  ; 3*i  ; 3*i  ];
            col(((i0+1):(i0+9))') = [3*i-2; 3*i-1; 3*i  ; 3*i-2; 3*i-1; 3*i  ; 3*i-2; 3*i-1; 3*i  ];
            
            M_block = 1/(config.method.sigma*delta_t)*factor* ((omega')*omega) + 1/h * eye(3,3);
            val(((i0+1):(i0+9))') = reshape(M_block',9,1);
            
            % S_iip
            ip = Gamma.neigh(i,2); 
            row(((i0+10):(i0+12))') = [3*i-2;  3*i-1;  3*i ];
            col(((i0+10):(i0+12))') = [3*ip-2; 3*ip-1; 3*ip];
            val(((i0+10):(i0+12))') = -1/h * [1; 1; 1];

            
            i0 = i0 + 12;
            
        else
            h = norm(Gamma.X_old(i,:)-Gamma.X_old(Gamma.neigh(i,1),:));
            factor = h/2;
            
            % S_ii
            row(((i0+1):(i0+9))') = [3*i-2; 3*i-2; 3*i-2; 3*i-1; 3*i-1; 3*i-1; 3*i  ; 3*i  ; 3*i  ];
            col(((i0+1):(i0+9))') = [3*i-2; 3*i-1; 3*i  ; 3*i-2; 3*i-1; 3*i  ; 3*i-2; 3*i-1; 3*i  ];
            
            M_block = 1/(config.method.sigma*delta_t)*factor* ((omega')*omega) + 1/h * eye(3,3);
            val(((i0+1):(i0+9))') = reshape(M_block',9,1);
            
            % S_iim
            im = Gamma.neigh(i,1);
            
            row(((i0+10):(i0+12))') = [3*i-2;  3*i-1;  3*i ];
            col(((i0+10):(i0+12))') = [3*im-2; 3*im-1; 3*im];
            val(((i0+10):(i0+12))') = -1/h * [1; 1; 1];
        
            
            i0 = i0 + 12; 
        end
    end
       
end

%% Reduce all vectors to necessary length
row_P_geo = row_P_geo(1:ig,1); 
col_P_geo = col_P_geo(1:ig,1); 
val_P_geo = val_P_geo(1:ig,1); 

row_P_triple = row_P_triple(1:it,1); 
col_P_triple = col_P_triple(1:it,1); 
val_P_triple = val_P_triple(1:it,1); 

row = row(1:i0,1);
col = col(1:i0,1);
val = val(1:i0,1); 


P_geo    = sparse(row_P_geo,col_P_geo,val_P_geo,3*J,3*J); 
P_triple = sparse(row_P_triple,col_P_triple,val_P_triple,3*J,3*J); 
S = sparse(row,col,val,3*J,3*J); 

P = P_triple * P_geo;
S = P*S*(P');
end


