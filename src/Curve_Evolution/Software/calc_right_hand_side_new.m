function b_out = calc_right_hand_side_new(Gamma,Omega,Surface,config,Image)

J = size(Gamma.X,1);
dim = size(Gamma.X,2);
b = zeros(J,dim);
F = zeros(J,1);

PX=Gamma.X_old;

if(size(size(Image.data),2)==3)
    I0 = double([Image.data(:,:,1), Image.data(:,:,2), Image.data(:,:,3)]); 
else
    I0 = double(Image.data); 
end

for j=1:size(Gamma.index_info,1)
    for jj=1:Gamma.nr_curves(j,1)
        i = Gamma.index_info(j,jj);
        ii = 1;

        while(i~=Gamma.index_info(j,jj) || ii==1)

            if(dim==2)
                % Faster with C-code
                F(i,1) = funcf(Gamma.X_old(i,:),Omega.coeffs,config.method.lambda,Gamma.orient,I0,[Image.sizes,size(Gamma.index_info,1),size(Omega.coeffs,1),config.method.color,j]);
            else
                % Fcn using closest_simp and Image.data
                F(i,1) = func_f(i,Gamma,Image,Omega,config,j);
            end

            
            % no external force at triple points because they belong to
            % more than one curve
            if(Gamma.neigh(i,1)==-2 || Gamma.neigh(i,2)==-3)
                F(i,1)=0;
            end
            
            
            if(Gamma.neigh(i,1)>0 && Gamma.neigh(i,2)>0)  % no boundary index i
                h1 = norm(Gamma.X_old(i,:)-Gamma.X_old(Gamma.neigh(i,1),:));
                h2 = norm(Gamma.X_old(Gamma.neigh(i,2),:)-Gamma.X_old(i,:));
                
                b(i,:) = 1/h1*PX(Gamma.neigh(i,1),:)+ 1/h2*PX(Gamma.neigh(i,2),:) -(1/h1+1/h2)*PX(i,:)  + 1/config.method.sigma*Gamma.nu(i,:)*(h1+h2)/2*F(i,1);
                
                
                ii = ii+1;
                i  = Gamma.neigh(i,2);
            else
                if(Gamma.neigh(i,1)<0)
                    h = norm(Gamma.X_old(Gamma.neigh(i,2),:)-Gamma.X_old(i,:));
                    
                    
                    b(i,:)= 1/h*PX(Gamma.neigh(i,2),:) - 1/h*PX(i,:) + 1/config.method.sigma*Gamma.nu(i,:)*h/2*F(i,1);
                    
                    
                    ii=ii+1;
                    i = Gamma.neigh(i,2);
                else  % Gamma.neigh(i,2)<0
                    h = norm(Gamma.X_old(i,:)-Gamma.X_old(Gamma.neigh(i,1),:));
                                       
                    b(i,:) = 1/h*PX(Gamma.neigh(i,1),:) - 1/h*PX(i,:) + 1/config.method.sigma*Gamma.nu(i,:)*h/2*F(i,1);
                                        
                    ii=ii+1;
                    i = Gamma.index_info(j,jj);  % end loop if boundary index is reached
                end
            end
        end
        
        
        
        
    end
end

if(dim==2)
    for k=1:size(Gamma.Lambda,1)
        nodes = Gamma.Lambda(k,1:3);
        b(nodes(1),:) = b(nodes(1),:) + b(nodes(2),:) + b(nodes(3),:);
        b(nodes(2),:) = [0,0];
        b(nodes(3),:) = [0,0];
    end
end


% Reshape b
b_out = reshape(b', dim*J,1); 


end


