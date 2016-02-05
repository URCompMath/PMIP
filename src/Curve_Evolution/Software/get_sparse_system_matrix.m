function S = get_sparse_system_matrix(Gamma,config,delta_t)

J=size(Gamma.X,1); 

% Sparse Matrix S
row = zeros(6*J,1); 
col = zeros(6*J,1); 
val = zeros(6*J,1); 

i0 = 0; 


N_triple = size(Gamma.Lambda,1); 
ind = [(1:J)',zeros(J,1)]; 
for i=1:N_triple
    ind(Gamma.Lambda(i,1),2) = i;  % Storing the Lambda index
    for j=2:3
        ind(Gamma.Lambda(i,j),1) = Gamma.Lambda(i,1); 
    end
end

for i=1:J
    omega = Gamma.nu(i,:)';

    if(Gamma.neigh(i,1) > 0 && Gamma.neigh(i,2) > 0)
        h1 = norm(Gamma.X_old(i,:)-Gamma.X_old(Gamma.neigh(i,1),:));
        h2 = norm(Gamma.X_old(Gamma.neigh(i,2),:)-Gamma.X_old(i,:)); 
        factor = (h1+h2)/2;
              
        % S_ii
        row(((i0+1):(i0+4))') = [2*i-1; 2*i-1; 2*i; 2*i];
        col(((i0+1):(i0+4))') = [2*i-1; 2*i; 2*i-1; 2*i];
        val(((i0+1):(i0+4))') =  1/(config.method.sigma*delta_t)*factor*[omega(1)^2; omega(1)*omega(2); omega(2)*omega(1); omega(2)^2] + (1/h1+1/h2)*[1;0;0;1];
    
        % S_iim and S_iip
        im = ind(Gamma.neigh(i,1),1); 
        ip = ind(Gamma.neigh(i,2),1); 
        row(((i0+5):(i0+8))') = [2*i-1; 2*i; 2*i-1; 2*i];
        col(((i0+5):(i0+8))') = [2*im-1;2*im;2*ip-1;2*ip];
        val(((i0+5):(i0+8))') = [-1/h1;-1/h1;-1/h2;-1/h2]; 
        
        i0 = i0+8;
    else
        if(Gamma.neigh(i,1)<0)
            if(i==ind(i,1))
                h=norm(Gamma.X_old(Gamma.neigh(i,2),:)-Gamma.X_old(i,:));
                factor = h/2;
                
                % S_ii
                row(((i0+1):(i0+4))') = [2*i-1; 2*i-1; 2*i; 2*i];
                col(((i0+1):(i0+4))') = [2*i-1; 2*i; 2*i-1; 2*i];
                val(((i0+1):(i0+4))') =  1/(config.method.sigma*delta_t)*factor*[omega(1)^2; omega(1)*omega(2); omega(2)*omega(1); omega(2)^2] + 1/h*[1;0;0;1];
                
                % S_iip
                ip = ind(Gamma.neigh(i,2),1);
                row(((i0+5):(i0+6))')=[2*i-1; 2*i];
                col(((i0+5):(i0+6))')=[2*ip-1;2*ip];
                val(((i0+5):(i0+6))')=-1/h*[1;1];
                
                
                % If representant of triple junctions 
                if(ind(i,2) > 0)
                    nodes = Gamma.Lambda(ind(i,2),2:3); 
                    for j=1:2
                        omega1 = Gamma.nu(nodes(j),:)'; 
                        if(Gamma.neigh(nodes(j),1) > 0)
                            h = norm(Gamma.X_old(nodes(j),:)-Gamma.X_old(Gamma.neigh(nodes(j),1),:)); 
                            ip = ind(Gamma.neigh(nodes(j),1),1);
                        else
                            h = norm(Gamma.X_old(nodes(j),:)-Gamma.X_old(Gamma.neigh(nodes(j),2),:)); 
                            ip = ind(Gamma.neigh(nodes(j),2),1);

                        end
                        
                        val(((i0+1):(i0+4))') = val(((i0+1):(i0+4))') + 1/(config.method.sigma*delta_t)*h/2*[omega1(1)^2; omega1(1)*omega1(2); omega1(2)*omega1(1); omega1(2)^2] + 1/h*[1;0;0;1];        
                        
                        row(((i0+6+2*j-1):(i0+6+2*j))')=[2*i-1; 2*i];
                        col(((i0+6+2*j-1):(i0+6+2*j))')=[2*ip-1;2*ip];
                        val(((i0+6+2*j-1):(i0+6+2*j))')=-1/h*[1;1];

                    end
                    i0 = i0 + 4; 
                end
                i0 = i0 + 6;
            else
                
                row(((i0+1):(i0+4))') = [2*i-1; 2*i-1;        2*i; 2*i];
                col(((i0+1):(i0+4))') = [2*i-1; 2*ind(i,1)-1; 2*i; 2*ind(i,1)];
                val(((i0+1):(i0+4))') = [1;-1;1;-1];
                i0 = i0 + 4; 
            end
                
            
        else
            if(i==ind(i,1))
                h = norm(Gamma.X_old(i,:)-Gamma.X_old(Gamma.neigh(i,1),:));
                factor = h/2;
                
                % S_ii
                row(((i0+1):(i0+4))') = [2*i-1; 2*i-1; 2*i; 2*i];
                col(((i0+1):(i0+4))') = [2*i-1; 2*i; 2*i-1; 2*i];
                val(((i0+1):(i0+4))') =  1/(config.method.sigma*delta_t)*factor*[omega(1)^2; omega(1)*omega(2); omega(2)*omega(1); omega(2)^2] + 1/h*[1;0;0;1];
                
                % S_iip
                im = ind(Gamma.neigh(i,1),1);
                row(((i0+5):(i0+6))')=[2*i-1; 2*i];
                col(((i0+5):(i0+6))')=[2*im-1;2*im];
                val(((i0+5):(i0+6))')=-1/h*[1;1];
                
                
                if(ind(i,2) > 0)
                    nodes = Gamma.Lambda(ind(i,2),2:3); 
                    for j=1:2
                        omega1 = Gamma.nu(nodes(j),:)'; 
                        if(Gamma.neigh(nodes(j),1) > 0)
                            h = norm(Gamma.X_old(nodes(j),:)-Gamma.X_old(Gamma.neigh(nodes(j),1),:)); 
                            ip = ind(Gamma.neigh(nodes(j),1),1);
                        else
                            h = norm(Gamma.X_old(nodes(j),:)-Gamma.X_old(Gamma.neigh(nodes(j),2),:)); 
                            ip = ind(Gamma.neigh(nodes(j),2),1);

                        end
                        
                        val(((i0+1):(i0+4))') = val(((i0+1):(i0+4))') + 1/(config.method.sigma*delta_t)*h/2*[omega1(1)^2; omega1(1)*omega1(2); omega1(2)*omega1(1); omega1(2)^2] + 1/h*[1;0;0;1];        
                        row(((i0+6+2*j-1):(i0+6+2*j))')=[2*i-1; 2*i];
                        col(((i0+6+2*j-1):(i0+6+2*j))')=[2*ip-1;2*ip];
                        val(((i0+6+2*j-1):(i0+6+2*j))')=-1/h*[1;1];

                    end
                    i0 = i0 + 4; 
                end
                
                i0 = i0 + 6;
            else
                
                row(((i0+1):(i0+4))') = [2*i-1; 2*i-1;        2*i; 2*i];
                col(((i0+1):(i0+4))') = [2*i-1; 2*ind(i,1)-1; 2*i; 2*ind(i,1)];
                val(((i0+1):(i0+4))') = [1;-1;1;-1];
                i0 = i0 + 4; 
            end
        end
    end
end



S = sparse(row,col,val,2*J,2*J);

end


           