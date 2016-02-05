function [mark,simplices_new] = mark_simplices_for_deletion(X,simplices,p0,n0,nr_surfaces,config)

N_simplices = size(simplices,1);

% Initialize output
mark = zeros(N_simplices,1);
simplices_new = simplices;

n0=n0/norm(n0);
tol_distn = config.topchange.a_factor*config.topchange.a;
if(config.image.flag == 9)  % set tol_distn smaller
    tol_distn = 6;
end

set1 = zeros(100,1); 
set2 = zeros(100,1); 
n1 = 0; 
n2 = 0; 

for i=1:N_simplices
    q0 = X(simplices{i,1}.nodes(1),:)';
    q1 = X(simplices{i,1}.nodes(2),:)';
    q2 = X(simplices{i,1}.nodes(3),:)';
    
    distn_1 = abs(n0'*(q0-p0));
    distn_2 = abs(n0'*(q1-p0));
    distn_3 = abs(n0'*(q2-p0));
    
    distn = min([distn_1, distn_2, distn_3]);
    
    
    if(distn < tol_distn * 4) 
        % Consider sign lambda_i to decide if simplex shoud be deleted
        lambda_1 = n0'*(q0-p0)/norm(q0-p0);   % = cos(alpha_1)
        lambda_2 = n0'*(q1-p0)/norm(q1-p0);   % = cos(alpha_2)
        lambda_3 = n0'*(q2-p0)/norm(q2-p0);   % = cos(alpha_3)
        
        s1 = sign(lambda_1);
        s2 = sign(lambda_2);
        s3 = sign(lambda_3);
        
        tol = cos((90-1)*pi/180); % alpha < 1° --> approx. on intersection plane

        lambda = min([lambda_1,lambda_2,lambda_3]);
        
        % Get index
        if(s1==s2 && s2==s3 && abs(lambda)> tol && distn > tol_distn)
            % No intersection of simplex with separation plane
            % Change sub-index in case of splitting
            if(s1 < 0)
                simplices_new{i,1}.index(2) = simplices{i,1}.index(2);
                n1 = n1+1; 
                set1(n1,1) = i; 
            else
                simplices_new{i,1}.index(2) = nr_surfaces + 1;
                n2 = n2 + 1; 
                set2(n2,1) = i; 
            end
        else
            % Mark simplex i for deletion
            mark(i)=1;
        end
    else
        simplices_new{i,1}.index(2) = -1; 
    end
end

%% Assign remaining simplices to a surface index by heritage
set1 = set1(1:n1,1); 
set2 = set2(1:n2,1); 

index2 = simplices_new{set1(1,1),1}.index(2);
while(n1>0)
    set1new = zeros(n1,1); 
    n1new = 0; 
    
    for k=1:n1
        i = set1(k,1); 
        
        for j=1:3
            ineigh = simplices_new{i,1}.neigh(j); 
            
            if(simplices_new{ineigh,1}.index(2)==-1)
                simplices_new{ineigh,1}.index(2) = index2;
                n1new = n1new +1; 
                set1new(n1new,1) = ineigh;
            end
        end
    end
    set1new = set1new(1:n1new,1); 
    
    set1 = set1new;
    n1 = size(set1,1); 
end
    
index22 = simplices_new{set2(1,1),1}.index(2);
while(n2>0)
    set2new = zeros(n2,1); 
    n2new = 0; 
    
    for k=1:n2
        i = set2(k,1); 
        
        for j=1:3
            ineigh = simplices_new{i,1}.neigh(j); 
            
            if(simplices_new{ineigh,1}.index(2)==-1)
                simplices_new{ineigh,1}.index(2) = index22;
                n2new = n2new +1; 
                set2new(n2new,1) = ineigh;
            end
        end
    end
    set2new = set2new(1:n2new,1); 
    
    set2 = set2new;
    n2 = size(set2,1); 
end    
    
    
    
end