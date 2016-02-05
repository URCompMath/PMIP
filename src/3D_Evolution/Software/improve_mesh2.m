function [Gamma_new,n_remaining] = improve_mesh2(Gamma,i_p0)
% Improve Gamma at i_p0 (to many edges end at i_p0)

%% 1. Search for simplices containing i_p0
list = zeros(100,2); % index simplex, local index node i_p0
n = 0; 

Nsimp = size(Gamma.simplices); 
for i=1:Nsimp
    for j=1:3
        if(Gamma.simplices{i,1}.nodes(j) == i_p0)
            n=n+1;
            list(n,:)=[i,j]; 
        end
    end
end
list = list(1:n,:);

%% 2. If n is odd, divide simplex with longest edge
ind_nr = [2,3; 3,1; 1,2]; 
if(2*ceil(n/2) ~= n) % n is odd
    % find largest opposite length
    i_longest = 0; 
    l_longest = -1;
    
    for i=1:n
        isimp = list(i,1); 
        e = Gamma.simplices{isimp,1}.edges(list(i,2)); 
        ix1 = Gamma.edges(e,1); 
        ix2 = Gamma.edges(e,2); 
        
        l = norm(Gamma.X(ix1,:)-Gamma.X(ix2,:)); 
        if(l > l_longest)
            l_longest = l;
            i_longest = i;
        end
    end
    
    % Divide simplex i_longest
    isimp = list(i_longest,1);
    kloc  = list(i_longest,2); % local index of i_p0
    
        
    i1 = Gamma.simplices{isimp,1}.nodes(ind_nr(kloc,1));
    i2 = Gamma.simplices{isimp,1}.nodes(ind_nr(kloc,2));
    e12 = Gamma.simplices{isimp,1}.edges(kloc) ;
    
    ineigh = Gamma.simplices{isimp,1}.neigh(kloc);
        
    found = 0; 
    j=1; 
    while(j<=3 && found==0)
        if(Gamma.simplices{ineigh,1}.edges(j)==e12)
            found = j; 
        end
        j=j+1;
    end
    kloc_neigh = found;
    i3 = Gamma.simplices{ineigh,1}.nodes(kloc_neigh);
    
    % a) One new node
    Xnew = (Gamma.X(i1,:) + Gamma.X(i2,:))/2;
    Gamma.X = [Gamma.X; Xnew]; 
    Nnodes = size(Gamma.X,1);
    
    % b) Three new edges
    e1 = [Nnodes, i3];
    e2 = [Nnodes, i2];
    e3 = [Nnodes, i_p0];
    
    Gamma.edges = [Gamma.edges; e1; e2; e3];   
    Nedges = size(Gamma.edges,1); % e1: Nedges-2, 
                                  % e2: Nedges-1,
                                  % e3: Nedges
                                  
    % c) Change edge e12
    Gamma.edges(e12,:) = [Nnodes,i1]; 
    
    
    % Two new simplices
    Nsimp = size(Gamma.simplices,1) + 2;
        
    
    % d) Change simplex isimp (nodes, edges,neigh)
    neigh_save1 = Gamma.simplices{isimp,1}.neigh(ind_nr(kloc,1));
    edges_save1 = Gamma.simplices{isimp,1}.edges(ind_nr(kloc,1));

    Gamma.simplices{isimp,1}.nodes(ind_nr(kloc,2)) = Nnodes;  % former i2, other nodes untouched
    Gamma.simplices{isimp,1}.edges(ind_nr(kloc,1)) = Nedges;  % e3 = Nedges, other edges untouched
    
    
    Gamma.simplices{isimp,1}.neigh(ind_nr(kloc,1)) = Nsimp;   % former neigh_save1, other neighbor relations untouched    
    
    
    % e) Change neigh info of former neighbor of isimp: neigh_save1
    found = 0; 
    j=1; 
    while(j<=3 && found==0)
        if(Gamma.simplices{neigh_save1,1}.neigh(j) == isimp)
            found = j;
        end
        j=j+1;
    end
    Gamma.simplices{neigh_save1,1}.neigh(found) = Nsimp;
    
    
    % f) Change simplex ineigh (nodes, edges,neigh)
    neigh_save2 = Gamma.simplices{ineigh,1}.neigh(ind_nr(kloc_neigh,2)); 
    edges_save2 = Gamma.simplices{ineigh,1}.edges(ind_nr(kloc_neigh,2));

    Gamma.simplices{ineigh,1}.nodes(ind_nr(kloc_neigh,1)) = Nnodes;   % former i2, other nodes untouched
    Gamma.simplices{ineigh,1}.edges(ind_nr(kloc_neigh,2)) = Nedges-2; % e1 = Nedges - 2, other edges untouched
    
    
    Gamma.simplices{ineigh,1}.neigh(ind_nr(kloc_neigh,2)) = Nsimp-1;  % former neigh_save2, other neighbor relations untouched
        
    
    % g) Change neigh info of former neighbor of ineigh: neigh_save2
    found = 0;
    j=1; 
    while(j<=3 && found==0)
        if(Gamma.simplices{neigh_save2,1}.neigh(j) == ineigh)
            found = j;
        end
        j=j+1;
    end
    Gamma.simplices{neigh_save2,1}.neigh(found) = Nsimp-1;
    
    
    % h) Create simplex Nsimp-1
    Gamma.simplices{Nsimp-1,1}.nodes = [i3, i2, Nnodes]; 
    Gamma.simplices{Nsimp-1,1}.edges = [Nedges-1, Nedges-2, edges_save2];  % e2, e1, edges_save2
    Gamma.simplices{Nsimp-1,1}.neigh = [Nsimp, ineigh, neigh_save2]; 
    Gamma.simplices{Nsimp-1,1}.index = Gamma.simplices{ineigh,1}.index;
        
    
    % i) Create simplex Nsimp
    Gamma.simplices{Nsimp,1}.nodes   = [i2,i_p0,Nnodes]; 
    Gamma.simplices{Nsimp,1}.edges   = [Nedges, Nedges-1, edges_save1]; % e3, e2, edges_save1
    Gamma.simplices{Nsimp,1}.neigh   = [isimp, Nsimp-1, neigh_save1]; 
    Gamma.simplices{Nsimp,1}.index   = Gamma.simplices{isimp,1}.index; 
    
    
    % j) Adapt list: add simplex Nsimp to list, i_p0 is stored in local
    % index 2 
    n = n+1;
    list(n,:) = [Nsimp, 2]; 
    

    
end


%% 3. Get ordered list 
list = [list, zeros(n,2)]; % i, k, i1, i2
for i=1:n
    i1 = Gamma.simplices{list(i,1),1}.nodes(ind_nr(list(i,2),1)); 
    i2 = Gamma.simplices{list(i,1),1}.nodes(ind_nr(list(i,2),2));
    list(i,3) = i1;
    list(i,4) = i2;
end
% list

list_new = zeros(n,4); 
list_new(1,:) = list(1,:); 

for i=2:n
    node_search = list_new(i-1,4); 
    index_found = find(list(:,3)==node_search);
    list_new(i,:) = list(index_found,:); 
end
if(size(list_new,1)~= n)
    fprintf('Error number of simplices in list_new %d is not equal to n=%d\n', size(list_new,1),n);
end
list = list_new;

%% 4. Create new nodes, edges, simplices
Nnodes = size(Gamma.X,1);
Nsimp  = size(Gamma.simplices,1);
Nedges = size(Gamma.edges,1);


% Create n/2 new nodes
for i=1:(n/2)
    i1 = list(2*i-1,3); 
    Gamma.X(Nnodes+i,:) = (Gamma.X(i1,:) + Gamma.X(i_p0,:))/2;
end

% Edges and simplices
for i=1:(n/2)
    % Change former edges
    isimp1 = list(2*i-1,1); 
    isimp2 = list(2*i,1); 
    kloc1 = list(2*i-1,2);  % local index of i_p0 in simplex isimp1
    kloc2 = list(2*i,2);    % local index of i_p0 in simplex isimp2
    
    e1 = Gamma.simplices{isimp1,1}.edges(ind_nr(kloc1,2)); 
    e2 = Gamma.simplices{isimp1,1}.edges(ind_nr(kloc1,1)); 
    e3 = Gamma.simplices{isimp2,1}.edges(ind_nr(kloc2,1)); 
    
    i1 = Gamma.simplices{isimp1,1}.nodes(ind_nr(kloc1,1)); 
    i2 = Gamma.simplices{isimp1,1}.nodes(ind_nr(kloc1,2)); 
    i3 = Gamma.simplices{isimp2,1}.nodes(ind_nr(kloc2,2)); 
    
    Gamma.edges(e1,:) = [i1, Nnodes+i]; 
    
    if(i<(n/2))
        ip = i+1;
    else
        ip = 1;
    end
    Gamma.edges(e2,:) = [Nnodes+i,Nnodes+ip]; 
    Gamma.edges(e3,:) = [i3, Nnodes+ip]; 
    
    % Create 3 new edges (per i)
    Gamma.edges(Nedges+3*i-2,:) = [i_p0, Nnodes+i]; 
    Gamma.edges(Nedges+3*i-1,:) = [i2, Nnodes+i]; 
    Gamma.edges(Nedges+3*i,:)   = [i2, Nnodes+ip]; 
    % set by next simplex: Gamma.edges(Nedges+3*(i-1)-2,:) = [i_p0,Nnodes+ip]
    
    % Change simplex isimp1
    Gamma.simplices{isimp1,1}.nodes(kloc1) = Nnodes+i;                   % former i_p0, other nodes untouched
    Gamma.simplices{isimp1,1}.edges(ind_nr(kloc1,1)) = Nedges + 3*i-1;   % former e2, other edges untouched
    Gamma.simplices{isimp1,1}.neigh(ind_nr(kloc1,1)) = Nsimp + 2*i-1;    % former isimp2, other neighbors untouched
    
    % Change simplex isimp2
    Gamma.simplices{isimp2,1}.nodes(kloc2) = Nnodes + ip;                % former i_p0, other nodes untouched
    Gamma.simplices{isimp2,1}.edges(ind_nr(kloc2,2)) = Nedges + 3*i;     % former e2, other edges untouched
    Gamma.simplices{isimp2,1}.neigh(ind_nr(kloc2,2)) = Nsimp + 2*i-1;    % former isimp1, other neighbors untouched
    
    % Create two new simplices
    if(i>=2)
        im = i-1;
    else
        im = n/2;
    end
    
    Gamma.simplices{Nsimp+2*i-1,1}.nodes = [Nnodes+ip, Nnodes+i, i2];
    Gamma.simplices{Nsimp+2*i-1,1}.edges = [Nedges+3*i-1, Nedges+3*i, e2]; 
    Gamma.simplices{Nsimp+2*i-1,1}.neigh = [isimp1, isimp2, Nsimp+2*i]; 
    Gamma.simplices{Nsimp+2*i-1,1}.index = Gamma.simplices{isimp1,1}.index; 
    
    Gamma.simplices{Nsimp+2*i,  1}.nodes = [Nnodes+i, Nnodes+ip, i_p0]; 
    Gamma.simplices{Nsimp+2*i,  1}.edges = [Nedges+3*ip-2, Nedges+3*i-2, e2]; 
    Gamma.simplices{Nsimp+2*i,  1}.neigh = [Nsimp+2*ip, Nsimp+2*im, Nsimp+2*i-1];
    Gamma.simplices{Nsimp+2*i,  1}.index = Gamma.simplices{isimp1,1}.index; 
   
    
end
    


%% Output
Gamma_new = Gamma; 
n_remaining = n/2; 
end
