function [index_klm,A,Xindex,omega] = get_indices_top_change_four_cases(Gamma,config,a,sizes,cc,Ndetect)
%% Preparations
X = Gamma.X;
simplices = Gamma.simplices;

size_x = sizes(1);
size_y = sizes(2);
size_z = sizes(3);

% Auxiliary matrix A for splitting detection 
Nx = 2*ceil(size_x/a);
Ny = 2*ceil(size_y/a);
Nz = 2*ceil(size_z/a);
A  = cell(Nx,Ny,Nz); 

max_nodes = 0;

index_klm = zeros(size(X,1),3); 
index_klm_nr = 0; 

% Preparation get surface index for each node (two phase, one main surface)
Xindex = zeros(size(X,1),1); 
for i=1:size(simplices,1)
    nodes = simplices{i,1}.nodes;
    index_sub = simplices{i,1}.index(2); 
    
    for k=1:3
        Xindex(nodes(k),1) = index_sub;
    end
end

% Preparation discrete normal vectors at nodes
omega = zeros(size(X,1),3);

for i=1:size(simplices,1)
    nodes = simplices{i,1}.nodes;
    for k=1:3
        % Update omega
        omega(nodes(k),:) = omega(nodes(k),:) + Gamma.area_sigma(i)*Gamma.nu_sigma(i,:);
    end
end
for k=1:size(X,1)
    omega(k,:) = omega(k,:) / norm(omega(k,:));
end


%% Main loop
for i=1:size(X,1)
    % Part I
    help = X(i,:)-cc;
    k = ceil(help(1)/a);
    l = ceil(help(2)/a);
    m = ceil(help(3)/a);
    

    
    if(k<=0)
        k=1;
    else
        if(k> Nx)
            k=Nx;
        end
    end
    if(l<=0)
        l=1;
    else
        if(l> Ny)
            l=Ny;
        end
    end
    if(m<=0)
        m=1;
    else
        if(m>Nz)
            m=Nz;
        end
    end
    
    
    if(isempty(A{k,l,m})==1)
        A{k,l,m}=struct('nr', 1, 'nodes', zeros(100,1), 'stored', 0, 'index_sub', zeros(100,1));
        A{k,l,m}.nodes(1,1)=i;
        A{k,l,m}.index_sub(1,1)= Xindex(i,1); 
    else 
        A{k,l,m}.nr = A{k,l,m}.nr + 1; 
        
        if(A{k,l,m}.nr > max_nodes) 
            max_nodes = A{k,l,m}.nr ;
        end
        A{k,l,m}.nodes(A{k,l,m}.nr,1)=i;
        A{k,l,m}.index_sub(A{k,l,m}.nr,1)=Xindex(i,1); 
        
        % Look for topological changes
        
        % 1. Look for different sub-indices
        isub = A{k,l,m}.index_sub(1,1); 
        if(Xindex(i,1) ~= isub)
            if(A{k,l,m}.stored == 0)
                A{k,l,m}.stored = 1; 
                index_klm_nr = index_klm_nr + 1; 
                index_klm(index_klm_nr,:) = [k,l,m]; 
                fprintf('Two different surface sub-indices --> top change detected\n'); 
            else
                A{k,l,m}.stored = 1; % mark that two different sub-indices occur
            end
        else
            % 2. Look for nr > Ndetect (--> probably splitting)
            if(A{k,l,m}.nr>=Ndetect)
                if(A{k,l,m}.stored == 0)
                    A{k,l,m}.stored = 2;
                    index_klm_nr = index_klm_nr + 1;
                    index_klm(index_klm_nr,:)=[k,l,m];
                    fprintf('Number of nodes %d > %d = Ndetect in a cube --> probable top change detected\n', A{k,l,m}.nr, Ndetect); 
                end
            else
                % 3. Consider normal vectors, if nodes belong to the same
                % surface but the local surface parts have opposite normal
                % vectors a merging takes place
                found = 0; 
                j=1; 
                while(j<= A{k,l,m}.nr-1 && found == 0)
                    omega1 = omega(A{k,l,m}.nodes(j,1),:); 
                    omega2 = omega(i,:);
                    alpha = acos(omega1*omega2');
                    if(alpha > config.topchange.merge_angle)
                        found = 1;
                        A{k,l,m}.stored = 3;
                        index_klm_nr = index_klm_nr + 1;
                        index_klm(index_klm_nr,:) = [k,l,m];
                        
                        fprintf('Two opposite normal vectors --> top change detected\n');
                        fprintf('alpha = %3.3f omega1 = (%2.3f %2.3f %2.3f)  omega2 = (%2.3f %2.3f %2.3f)\n', alpha*180/pi, omega1(1),omega1(2),omega1(3),omega2(1),omega2(2),omega2(3)); 
                        
                    end
                        
                    j=j+1;
                end
            end
        
        end
            
            
    end
    
end

% Reduce matrix rows
index_klm = index_klm(1:index_klm_nr,:);



%% Printing nodes/result of get_indices_top_change
if(index_klm_nr > 0)
    fprintf('Squares where a top change likely occurs:\n'); 
end

for i=1:index_klm_nr
    k = index_klm(i,1); 
    l = index_klm(i,2); 
    m = index_klm(i,3); 
    fprintf('%d/%d\nSquare %d %d %d, Stored index = %d,\n', i, index_klm_nr, k,l,m, A{k,l,m}.stored); 
    fprintf('Nodes:\n'); 
    for j=1:A{k,l,m}.nr
        k_node  = A{k,l,m}.nodes(j,1); 
        fprintf('node nr = %d  nodes = (%2.3f %2.3f %2.3f)\n', k_node, X(k_node,1), X(k_node,2), X(k_node,3)); 
    end
    fprintf('\n\n'); 
end


end
