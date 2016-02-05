function [X_new,simplices_new,edges_new,i_p0] = generate_new_simplices2(X,simplices,edges,p0)
% nodes ordered clockwise per simplex!


% Update X
X_new = [X; p0'; p0']; 
Nnodes = size(X_new,1);
i_p0 = [Nnodes-1, Nnodes];

% Update simplices and edges
Nsimp = size(simplices,1); 

info = zeros(500,5); % i, k, nodes(ind_nr(k,1)), nodes(ind_nr(k,1)), 0
ninfo = 0; 
 

% Get simplices with free edges
ind_nr = [2,3; 3,1; 1,2]; 
for i=1:Nsimp    
    for k=1:3
        if(simplices{i,1}.neigh(k)==-1)
                ninfo = ninfo + 1;
                info(ninfo,:) = [i,k, simplices{i,1}.nodes(ind_nr(k,1)), simplices{i,1}.nodes(ind_nr(k,2)), 0]; 
        end
    end
end
info = info(1:ninfo,:); 

% Get info1 and info2
info1 = zeros(ninfo,5); 

info1(1,:) = info(1,:); 
info1(1,5) = Nnodes-1;
info(1,5) = 1; 

i=2;
node_start =info1(1,3);
node = 0;
while(i==2 || node~= node_start)
    node_search = info1(i-1,4); 
    node_index = find(info(:,3)==node_search);
    
    info1(i,:) = info(node_index,:); 
    info1(i,5) = Nnodes-1; 
    
    info(node_index,5) = 1;
    
    node = info1(i,4); 
    i=i+1;
end
info1 = info1(1:(i-1),:);
ninfo1 = size(info1,1); 

% info
remain_ind = find(info(:,5)==0);

i=remain_ind(2); 
info2 = zeros(ninfo,5); 

info2(1,:) = info(i,:); 
info2(1,5) = Nnodes;
info(i,5) = 2; 

i=2;
node_start = info2(1,3); 
node = 0; 
while(i==2 || node ~= node_start)
    node_search = info2(i-1,4); 
    node_index = find(info(:,3)==node_search); 
    
    info2(i,:) = info(node_index,:); 
    info2(i,5) = Nnodes;
    
    info(node_index,5) = 2; 
    
    node = info2(i,4); 
    i=i+1; 
end
info2 = info2(1:(i-1),:);
ninfo2 = size(info2,1); 


%% Part 1
Nsimp = size(simplices,1); 
Nedges = size(edges,1); 

% Create new edges
for i=1:ninfo1
    edges(Nedges+i,:) = [info1(i,3), info1(i,5)];
end


% Create new simplices
for i=1:ninfo1
    simplices{Nsimp+i,1}.nodes = [info1(i,4), info1(i,3), info1(i,5)]; 
    
    % Simplex in "plus" direction, neigh opposite of info1(i,3)
    if(i<ninfo1)
        ip = Nsimp+i+1;
    else
        ip = Nsimp+1;
    end
    % Simplex in "minus" direction, neigh opposite of info1(i,4)
    if(i>1)
        im = Nsimp+i-1;
    else
        im = Nsimp+ninfo1; 
    end
    
    % Edge between Nedges+i and ip
    if(i<ninfo1)
        ep = Nedges+i+1;
    else
        ep = Nedges+1;
    end
    
    % Edge between Nedges+i and im
    em = Nedges+i;
    
    % Modify neighbor information of info1(i,1) at former free edge
    simplices{info1(i,1),1}.neigh(info1(i,2)) = Nsimp+i; 
   
    
    % Create simplex "Nsimp+i"
    simplices{Nsimp+i,1}.neigh = [im, ip, info1(i,1)]; 
    simplices{Nsimp+i,1}.edges = [em, ep, simplices{info1(i,1),1}.edges(info1(i,2))];    
    simplices{Nsimp+i,1}.index = simplices{info1(i,1),1}.index;    
end


%% Part 2
Nsimp = size(simplices,1); 
Nedges = size(edges,1); 

% Create new edges
for i=1:ninfo2
    edges(Nedges+i,:) = [info2(i,3), info2(i,5)];
end


% Create new simplices
for i=1:ninfo2
    simplices{Nsimp+i,1}.nodes = [info2(i,4), info2(i,3), info2(i,5)]; 
    
    % Simplex in "plus" direction, neigh opposite of info2(i,3)
    if(i<ninfo2)
        ip = Nsimp+i+1;
    else
        ip = Nsimp+1;
    end
    % Simplex in "minus" direction, neigh opposite of info2(i,4)
    if(i>1)
        im = Nsimp+i-1;
    else
        im = Nsimp+ninfo2; 
    end
    
    % Edge between Nedges+i and ip
    if(i<ninfo2)
        ep = Nedges+i+1;
    else
        ep = Nedges+1;
    end
    
    % Edge between Nedges+i and im
    em = Nedges+i;
    
    % Modify neighbor information of info2(i,1) at former free edge
    simplices{info2(i,1),1}.neigh(info2(i,2)) = Nsimp+i; 
   
    
    % Create simplex "Nsimp+i"
    simplices{Nsimp+i,1}.neigh = [im, ip, info2(i,1)]; 
    simplices{Nsimp+i,1}.edges = [em, ep, simplices{info2(i,1),1}.edges(info2(i,2))];    
    simplices{Nsimp+i,1}.index = simplices{info2(i,1),1}.index; 
end


%% Output
edges_new = edges;
simplices_new = simplices; 

       
end