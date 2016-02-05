function Gamma_new = close_open_hole(Gamma,simp_free,free_nodes)
% Input:
% Gamma                 struct, structure of current surface
% simp_free             N x 3 double, list of free simplices, 
%                       1st col: simplex index, 2nd col: local edge index
%                       3rd col: group number (nominal = 1)
% free_nodes            N x 2 double, list of free nodes corresponding to
%                       the free simplices (nodes at the free edge)

% Output
% Gamma_new             struct, updated surface, hole closed

% Get ordered list of free simplices and free_nodes
N = size(simp_free,1); 
i = 1; 
j = 1; 
search = free_nodes(:,1); 

% Create auxiliary list
simp_free_ordered = zeros(N,5); 

while(j<=N)
    % Copy from line i in simp_free to line j in simp_free_ordered
    simp_free_ordered(j,1:3) = simp_free(i,:); 
    % Copy free nodes to column 4:5
    simp_free_ordered(j,4:5) = free_nodes(i,:); 

    % Go to neighboring free edge,  update i
    node_search = free_nodes(i,2); 
    i = find(search == node_search); 
    
    % Increase j
    j=j+1;
end


% Update free_nodes to ordered list
free_nodes = simp_free_ordered(:,4:5); 
% Update simp_free
simp_free = simp_free_ordered(:,1:3); 



%% Change Gamma (add N new simplices, connect free nodes to the mid point)
Nsimp = size(Gamma.simplices,1); 
Nnodes = size(Gamma.X,1); 
Nedges = size(Gamma.edges,1); 

% Increase number of simplices, Nnodes, Nedges
Nsimp = Nsimp + N; 
Nnodes = Nnodes + 1; 
Nedges = Nedges + N; 

% Create new node in the middle of the hole (average of free nodes) and 
% update Gamma.X
x0 = zeros(1,3); 
for i=1:N
    x0 = x0 + Gamma.X(free_nodes(i,1),:);
end
x0 = x0/N; 
Gamma.X(Nnodes,:) = x0; 

% Get new edges
for i=1:N
    Gamma.edges(Nedges - N + i,:) = [free_nodes(i,1), Nnodes]; 
end

% Generate new simplices
% i=1
isimp = Nsimp - N + 1;
iparent = simp_free(1,1); 
iedge = Gamma.simplices{iparent,1}.edges(simp_free(1,2)); 

Gamma.simplices{isimp,1}.nodes = [free_nodes(1,2), free_nodes(1,1), Nnodes]; 
Gamma.simplices{isimp,1}.neigh = [Nsimp, Nsimp-N+2, iparent]; 
Gamma.simplices{isimp,1}.edges = [Nedges-N+1, Nedges-N+2, iedge]; 
Gamma.simplices{isimp,1}.index = Gamma.simplices{iparent,1}.index;

%i=2:(N-1)
for i=2:(N-1)
    isimp = Nsimp - N + i; 
    iparent = simp_free(i,1); 
    iedge = Gamma.simplices{iparent,1}.edges(simp_free(i,2)); 
    
    Gamma.simplices{isimp,1}.nodes = [free_nodes(i,2), free_nodes(i,1), Nnodes]; 
    Gamma.simplices{isimp,1}.neigh = [Nsimp-N+i-1, Nsimp-N+i+1, iparent]; 
    Gamma.simplices{isimp,1}.edges = [Nedges-N+i, Nedges-N+i+1, iedge]; 
    Gamma.simplices{isimp,1}.index = Gamma.simplices{iparent,1}.index;
end

%i=N
isimp = Nsimp; 
iparent = simp_free(N,1); 
iedge = Gamma.simplices{iparent,1}.edges(simp_free(N,2)); 

Gamma.simplices{isimp,1}.nodes = [free_nodes(N,2), free_nodes(N,1), Nnodes]; 
Gamma.simplices{isimp,1}.neigh = [Nsimp-1, Nsimp-N+1, iparent]; 
Gamma.simplices{isimp,1}.edges = [Nedges, Nedges-N+1, iedge]; 
Gamma.simplices{isimp,1}.index = Gamma.simplices{iparent,1}.index; 


%% Print information
fprintf('Nnodes  %d, Nsimp %d, Nedges %d\n\n\n', Nnodes,Nsimp,Nedges); 
fprintf('New edges:\n'); 
for i=1:N
    fprintf('edge %d: nodes %d %d\n', Nedges-N+i, Gamma.edges(Nedges-N+i,1), Gamma.edges(Nedges-N+i,2)); 
end
fprintf('\n\nNew simplices\n'); 
for i=1:N
    nodes = Gamma.simplices{Nsimp-N+i,1}.nodes;
    neigh = Gamma.simplices{Nsimp-N+i,1}.neigh;
    edges = Gamma.simplices{Nsimp-N+i,1}.edges;
    
    fprintf('simp %d, nodes = %d %d %d, neigh = %d %d %d, edges = %d %d %d\n', Nsimp-N+i, nodes(1),nodes(2),nodes(3),neigh(1),neigh(2),neigh(3),edges(1),edges(2),edges(3)); 
end
fprintf('\n'); 


%% Output
Gamma_new = Gamma;

end