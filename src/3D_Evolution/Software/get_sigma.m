function sigma = get_sigma(Gamma,sigma_t0,config)

% Initialize sigma
J = size(Gamma.X,1); 
sigma = sigma_t0 * ones(J,1); 


% Store for each node the weighted normal vectors w_j^m of the neighbors
w_info = cell(J,1); 
Nedges = size(Gamma.edges,1); 
for i=1:Nedges
    i1 = Gamma.edges(i,1); 
    i2 = Gamma.edges(i,2); 
    
    if(isempty(w_info{i1,1}))
        w_info{i1,1}.nr = 1;
        w_info{i1,1}.w_neighs = zeros(100,3); 
        w_info{i1,1}.w_neighs(1,:) = Gamma.omega(i2,:); 
    else
        w_info{i1,1}.nr = w_info{i1,1}.nr + 1;
        w_info{i1,1}.w_neighs(w_info{i1,1}.nr,:) = Gamma.omega(i2,:); 
    end
    if(isempty(w_info{i2,1}))
        w_info{i2,1}.nr = 1;
        w_info{i2,1}.w_neighs = zeros(100,3); 
        w_info{i2,1}.w_neighs(1,:) = Gamma.omega(i1,:); 
    else
        w_info{i2,1}.nr = w_info{i2,1}.nr + 1;
        w_info{i2,1}.w_neighs(w_info{i2,1}.nr,:) = Gamma.omega(i1,:); 
    end
end

% Set sigma(i,1) = 0 if large difference between w_k^m and w_j^m for a
% neighbor node j
alpha_max = config.mesh.alpha_max;
for i=1:J
    
    n = w_info{i,1}.nr;
    found = 0;
    j=1;
    wk = Gamma.omega(i,:); 
    wk = wk/norm(wk);
    while(j<= n && found == 0)
        
        wj = w_info{i,1}.w_neighs(j,:);
        wj = wj/norm(wj); 
        
        alpha = acos(wk*(wj')); 
        if(alpha > alpha_max)
            sigma(i,1) = 0;
            found = 1;
        end
        j=j+1;
    end
end
        
    
end