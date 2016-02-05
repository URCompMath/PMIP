function top_change_index = identify_top_change(nodes,Xindex,omega,X,config)
% Input
% nodes                   N x 1, index of nodes involved in the top change
% Xindex                  Nnodes x 1, sub-index of surface where nodes
%                         belong to
% omega                   Nnodes x 1, weighted normal vector
% X                       Nnodes x 3, coordinates of the nodes
% config                  structure, configuration parameters


% Output top_change_index, flag 
% 0:   no top change
% 1:   splitting
% 2:   merging
% 3:   increasing genus
% 4:   decreasing genus


top_change_index = 18; 

thr1 = config.topchange.thr(1); 
thr2 = config.topchange.thr(2); 
thr3 = config.topchange.thr(3); 

 

N = size(nodes,1); 

x = zeros(N,3); 
w = zeros(N,3); 
for i=1:N
    x(i,:) = X(nodes(i),:); 
    w(i,:) = omega(nodes(i),:); 
end

%% 1. Group finding
N1 = 1;
N2 = 0; 

group1 = zeros(N,1); 
group2 = zeros(N,1); 

group1(N1,1) = 1; % first normal vector forms one group

w1 = w(1,:); 
w2 = [0,0,0]; 

repeat = 1; 

while(repeat > 0 && repeat <= N)
    mark = zeros(N,1);
    mark(1,1) = 1;
    
    % First loop
    for i=2:N
        % Check if node belongs to group1
        wi = w(i,:);
        alpha = acos(w1*wi');
        if(alpha <= thr1)
            N1 = N1 + 1;
            group1(N1,1) = i;
            mark(i,1) = 1;
        else
            if(N2==0)
                % Create new groups if alpha >= thr2
                if(alpha >= thr2)
                    N2 = 1;
                    group2(N2,1) = i;
                    w2 = w(i,:);
                    mark(i,1) = 1;
                end
            else
                % Check if node belongs to group2
                beta = acos(w2*wi');
                if(beta <= thr1)
                    N2 = N2 + 1;
                    group2(N2,1) = i;
                    mark(i,1) = 1;
                end
            end
        end
    end
    
    % If group2 is empty and group1 contains nearly all nodes, no top change
    if(N2 == 0 && N1 >= 0.9*N)
        top_change_index = 0;
        return;
    end
    
    % If group1 has a small number of nodes, replace representative and repeat
    if(N1 < 0.1*N || N2 < 0.1*N)
        repeat = repeat + 1;
        N1 = 1;
        N2 = 0;
        gsave = group1(1,1);
        group1 = zeros(N,1);
        group2 = zeros(N,1);
        group1(N1,1) = gsave + 1;
        
        if(repeat <=N)
            w1 = w(repeat,:);
            w2 = [0,0,0];
        end
    else
        repeat = 0;
    end
   
end


% repeat
if(repeat == N+1)
    top_change_index = 0; 
    return;
end


% Next loops
Nsum = sum(mark); 
iter = 1; 
maxiter = 10; 

while(Nsum < N && iter <= maxiter)
    % average normal vector group 1
    w1 = [0,0,0]; 
    for j=1:N1
        w1 = w1 + w(group1(j,1),:); 
    end
    w1 = w1 / norm(w1);
    
    % average normal vector group 2
    w2 = [0,0,0]; 
    for j=1:N2
        w2 = w2 + w(group2(j,1),:); 
    end
    w2 = w2 / norm(w2); 
    
    % Loop over remaining indices
    for i=2:N
        if(mark(i,1)==0)
            wi = w(i,:); 
            alpha = acos(w1*wi'); 
            ind = 1; 
            if(N2 > 0)
                beta = acos(w2*wi'); 
                if(beta < alpha)
                    alpha = beta; 
                    ind = 2; 
                end
            end
            
            if(alpha < thr1)
                if(ind==1)
                    N1 = N1+1;
                    group1(N1,1) = i; 
                    mark(i,1) = 1; 
                else
                    N2 = N2+1; 
                    group2(N2,1) = i; 
                    mark(i,1) = 1; 
                end
            end
        end
    end
    
    Nsum = sum(mark); 
    iter = iter + 1; 
end
        
                    
%% 2. Compute angles of remaining nodes for third / forth group
% average normal vector group 1
w1 = [0,0,0];
for j=1:N1
    w1 = w1 + w(group1(j,1),:);
end
w1 = w1 / norm(w1);

% average normal vector group 2
w2 = [0,0,0];
for j=1:N2
    w2 = w2 + w(group2(j,1),:);
end
w2 = w2 / norm(w2);


n_large = 0; 
for i=2:N
    if(mark(i,1)==0)
        alpha = acos(w(i,:)*w1'); 
        beta  = acos(w(i,:)*w2'); 
        if(alpha >= thr3 && beta >= thr3)
            n_large = n_large + 1;
        end
    end
end



Nthr = N/4 ;
if(n_large > Nthr)
    % More than two groups of normal vectors: splitting or decreasing genus
    z = mean(x); % Mean of the nodes
    
    pointing_inside = zeros(N,1); 
    for i=1:N
        v = z - x(i,:); 
        v = v/norm(v); 
        wi = w(i,:); 
        alpha = acos(v*wi'); 
        if(alpha < pi/2)
            pointing_inside(i,1) = 1; 
        end
    end
    
    Ninside = sum(pointing_inside);
    ratio = Ninside / N;
    if(ratio > 0.7)
        fprintf('Splitting\n'); 
        top_change_index = 1;
        
    else
        if(ratio < 0.3)
            fprintf('Decreasing genus\n');
            top_change_index = 4;
            
        end
    end
    
else
    % Look at surface index
    ind1 = Xindex(nodes(1,1),1); 
    
    found = 0;  
    j=2;
    while(j<=N && found==0)
        indj = Xindex(nodes(j,1),1); 
        if(ind1 ~= indj)
            found = 1;
        end
        j=j+1;
    end
    if(found == 0)
        top_change_index = 3;
        fprintf('Increasing genus\n'); 
    else
        top_change_index = 2;
        fprintf('Merging\n'); 
    end
    
end
        





end
