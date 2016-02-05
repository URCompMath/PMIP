function [Gamma_new, N_nodes_new] = coarsen_curve(Gamma,node_i0,N_nodes)

J = size(Gamma.X,1); 
index = zeros(N_nodes,1);
info  = zeros(J,1);
count = 0;

i=node_i0;
j=1;
while(i~=node_i0 || j==1)
    if(j==2*round(j/2) && Gamma.neigh(i,2)>0) % j is even, do not delete boundary points or triple points
        count = count + 1;
        index(count)=i;
        info(i)=1;
    end
    
    i = Gamma.neigh(i,2);
    if(i<0)
        i=node_i0;
    end
    j=j+1;
end
index = index(1:count,1);
N_nodes_new = N_nodes - count;

% Change neighbor info
i=node_i0;
j=1;
while(i~=node_i0 || j==1)
    if(j~=2*round(j/2))  % j is uneven
        help = Gamma.neigh(i,2);
        if(Gamma.neigh(i,1)>0 && info(Gamma.neigh(i,1))==1)
            i1 = Gamma.neigh(i,1);
            Gamma.neigh(i,1) = Gamma.neigh(i1,1);
        end
        if(Gamma.neigh(i,2)>0 && info(Gamma.neigh(i,2))==1)
            i1 = Gamma.neigh(i,2);
            Gamma.neigh(i,2) = Gamma.neigh(i1,2);
        end
    end
    i = help;
    j = j+1;
    if(i<0)
        i=node_i0;
    end
end


% reduce index number of other nodes
for i=1:J
    j=Gamma.neigh(i,1);
    if(j>0)
        Gamma.neigh(i,1)=j-sum(info(1:j,1));
    end
    j=Gamma.neigh(i,2);
    if(j>0)
        Gamma.neigh(i,2)=j-sum(info(1:j,1));
    end
end

% reduce index number of leading nodes in index_info
Nmain = size(Gamma.index_info,1);
for k=1:Nmain
    for l=1:Gamma.nr_curves(k,1)
        j=Gamma.index_info(k,l);
        Gamma.index_info(k,l)=j-sum(info(1:j,1));
    end
end

% reduce index number of nodes stored in Lambda
N_triple = size(Gamma.Lambda,1);
for i=1:N_triple
    for ii=1:3
        j=Gamma.Lambda(i,ii);
        Gamma.Lambda(i,ii)=j-sum(info(1:j,1));
    end
end


% Delete rows in X and neigh 
Gamma.X(index,:)=[];
Gamma.neigh(index,:)=[];
% Delete rows in closest_simp (dim=3)
if(size(Gamma.X,2)==3)
    Gamma.closest_simp(index,:) = [];
end

Gamma_new = Gamma; 
end