function Gamma_new = get_length_info(Gamma,config,Image,Surface)

N_nodes_min = 2;

Nmain = size(Gamma.index_info,1);
Nsub  = size(Gamma.index_info,2);

% Initialize output
Gamma.length  = zeros(Nmain,Nsub);
Gamma.mark    = zeros(Nmain,Nsub);

iter_max = 5;

for k=1:Nmain
    for l=1:Gamma.nr_curves(k,1)
        repeat = 1;
        iter = 1;
        
        coarsening_allowed = 1;
        
        while(repeat > 0 && iter <= iter_max)

            repeat = 0;
            Gamma.length(k,l) = 0;
            
            i=1;
            j=Gamma.index_info(k,l);

            while(j~=Gamma.index_info(k,l) || i==1)

                if(Gamma.neigh(j,1)<0 && Gamma.neigh(j,2)<0) % Curve consists of only one node
                    Gamma.mark(k,l)=1;
                else
                    if(Gamma.neigh(j,2)>0)
                        h = norm(Gamma.X(Gamma.neigh(j,2),:)-Gamma.X(j,:));
                        Gamma.length(k,l)=Gamma.length(k,l) + h;
                        
                        if(h>config.length.Lmax) 
                            repeat = 1;
                            Gamma = refine_curve_loc(Gamma,Surface,j); % insert a point between j and neigh(j,2)
                            Gamma.X_old = Gamma.X;  % as number of nodes could have been changed
                            coarsening_allowed = 0; % do not allow coarsening in this step
                            j = Gamma.neigh(Gamma.neigh(j,2),2);  % jump 2 nodes forward
                            
                        else
                            j = Gamma.neigh(j,2);
                        end
                        
                        
                        i = i+1;
                    else
                        j = Gamma.index_info(k,l); % stop loop if boundary point is reached
                        
                    end
                end
            end
            
            iter = iter + 1;

        end
        
        N_nodes = i-1;
        
        tol = config.length.tol(1,1);
        if(Gamma.age(k,l)> config.length.age_min && Gamma.neigh(Gamma.index_info(k,l),1)>0) 
            tol = config.length.tol(1,2); % Bigger tolerance for "old" curves which are closed, new and open curves have a smaller tolerance! (Are not deleted that early)
        end
        
        
        if(Gamma.length(k,l)< tol || (N_nodes<=N_nodes_min && Gamma.neigh(Gamma.index_info(k,l),1)>=-1)) 
            Gamma.mark(k,l)=1;
            fprintf('Will delete Curve Nr. %d %d. Length = %2.2f, tol = %2.2f, N_nodes = %d, N_nodes_min = %d\n', k, l, Gamma.length(k,l), tol,N_nodes,N_nodes_min);
        else
            while(Gamma.length(k,l)/N_nodes > config.length.Lmax || N_nodes <= N_nodes_min*3)  % Refine curve if length/N_nodes is too big
                fprintf('Will refine Curve Nr. %d %d, length = %2.2f, n = %d\n', k,l,Gamma.length(k,l),N_nodes);
                [Gamma,N_nodes] = refine_curve(Gamma,Image,Surface,Gamma.index_info(k,l),N_nodes);
                Gamma.X_old = Gamma.X;  % as number of nodes could have been changed in refinement/coarsening
            end
            if(coarsening_allowed)
                while(Gamma.length(k,l)/N_nodes < config.length.Lmin && N_nodes>config.length.Nmin) % Coarsen curve but only if length/N_nodes is too small
                    fprintf('Will coarsen Curve Nr. %d %d, length = %2.2f, n = %d\n', k,l,Gamma.length(k,l),N_nodes);
                    [Gamma,N_nodes] = coarsen_curve(Gamma,Gamma.index_info(k,l),N_nodes);
                    Gamma.X_old = Gamma.X;  % as number of nodes could have been changed in refinement/coarsening
                end
            end
            
        end
        
        
    end
    
    
end
Gamma_new = Gamma;
end