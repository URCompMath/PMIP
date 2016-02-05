function points = compute_cross_sections(Gamma,val,dim_flag)
% Gamma: surface structure
% dim_flag: 1,2,3 --> fixed x,y,z
% val: value, for example: dim_flag = 1; val = 50, compute cross section of
% Gamma with plane {(val,:,:)}
%
% Output: points, Nnodes x 3

Nsimp = size(Gamma.simplices,1); 

points = zeros(Nsimp*2,3); 
npoints = 0; 

for i=1:Nsimp
    lower   = zeros(1,3);
    greater = zeros(1,3); 
    stop = 0;
    
    for j=1:3
        val_test = Gamma.X(Gamma.simplices{i,1}.nodes(j),dim_flag);
        if(val_test <= val)
            lower(j) = 1; 
        else
            if(val_test >= val)
                greater(j) = 1;
            end
        end
        if(lower(j)==1 && greater(j)==1)
            % Vertex on plane
            npoints = npoints +1; 
            points(npoints,:) = Gamma.X(Gamma.simplices{i,1}.nodes(j),:);
            stop = 1;
        end
    end
    
    if(stop == 0)
        if(sum(lower) >= 1 && sum(greater) >= 1)
            % Triangle is truncated
            % Since stop == 0, no vertex is on the separation plane -->
            % i.e. one vertex on one side, two vertices on the other side
            if(sum(lower)==1)
                k=find(lower==1); 
                k1 = mod(k,3)+1;
                k2 = mod(k1,3)+1;                
                j  = Gamma.simplices{i,1}.nodes(k); 
                j1 = Gamma.simplices{i,1}.nodes(k1); 
                j2 = Gamma.simplices{i,1}.nodes(k2); 
                
                % Intersection dim_flag == val with [X(j), X(j1)]
                % Line: lambda mapsto X(j,:) + lambda ( X(j,:) - X(j1,:)) 
                % Equation to solve: val = X(j,dim_flag) + lambda ( X(j1,dim_flag) - X(j,dim_flag)) 
                lambda = (val - Gamma.X(j,dim_flag))/(Gamma.X(j1,dim_flag)-Gamma.X(j,dim_flag)); 
                p = Gamma.X(j,:) + lambda * (Gamma.X(j1,:) - Gamma.X(j,:));
                npoints = npoints + 1; 
                points(npoints,:) = p; 
                
                % Line: lambda mapsto X(j,:) + lambda ( X(j,:) - X(j2,:)) 
                % Equation to solve: val = X(j,dim_flag) + lambda ( X(j2,dim_flag) - X(j,dim_flag)) 
                lambda = (val - Gamma.X(j,dim_flag))/(Gamma.X(j2,dim_flag)-Gamma.X(j,dim_flag)); 
                p = Gamma.X(j,:) + lambda * (Gamma.X(j2,:) - Gamma.X(j,:));
                npoints = npoints + 1; 
                points(npoints,:) = p; 
                
                
            else
                if(sum(greater == 1))
                    k = find(greater == 1); 
                    k1 = mod(k,3)+1;
                    k2 = mod(k1,3)+1; 
                    j  = Gamma.simplices{i,1}.nodes(k);
                    j1 = Gamma.simplices{i,1}.nodes(k1);
                    j2 = Gamma.simplices{i,1}.nodes(k2);
                    
                    % Intersection dim_flag == val with [X(j), X(j1)]
                    % Line: lambda mapsto X(j,:) + lambda ( X(j,:) - X(j1,:))
                    % Equation to solve: val = X(j,dim_flag) + lambda ( X(j1,dim_flag) - X(j,dim_flag))
                    lambda = (val - Gamma.X(j,dim_flag))/(Gamma.X(j1,dim_flag)-Gamma.X(j,dim_flag));
                    p = Gamma.X(j,:) + lambda * (Gamma.X(j1,:) - Gamma.X(j,:));
                    npoints = npoints + 1;
                    points(npoints,:) = p;
                    
                    % Line: lambda mapsto X(j,:) + lambda ( X(j,:) - X(j2,:))
                    % Equation to solve: val = X(j,dim_flag) + lambda ( X(j2,dim_flag) - X(j,dim_flag))
                    lambda = (val - Gamma.X(j,dim_flag))/(Gamma.X(j2,dim_flag)-Gamma.X(j,dim_flag));
                    p = Gamma.X(j,:) + lambda * (Gamma.X(j2,:) - Gamma.X(j,:));
                    npoints = npoints + 1;
                    points(npoints,:) = p;

                end
            end
                    
                    
        end
    end
end
points = points(1:npoints,:); 
end
            
    