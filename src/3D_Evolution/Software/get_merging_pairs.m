function [Gamma_new, assign_simp, warning] = get_merging_pairs(Gamma,simp_free)
% Input
% Gamma              struct, current surface structure
% simp_free          Nx3,  stores free simplices with one edge == -1 (free
%                    edge, cols: index simp, local index free edge, group
%                    index 

N = size(simp_free,1);

%% Get nr edges of the two groups and corner points of the edges
index_edges_corners1 = zeros(N,2); 
index_edges1 = zeros(N,1); 
index_simp1 = zeros(N,1); 
index_edges_corners2 = zeros(N,2); 
index_edges2 = zeros(N,1); 
index_simp2 = zeros(N,1); 


N1 = 0;
N2 = 0;
ind_nr = [2,3; 3,1; 1,2]; 
for i=1:N
    i0 = simp_free(i,3);
    
    isimp = simp_free(i,1);
    k0    = simp_free(i,2);
    e     = Gamma.simplices{isimp,1}.edges(k0);
    
    if(i0==1)
        N1 = N1 + 1;
        index_edges1(N1,1) = e;
        index_edges_corners1(N1,:) = [Gamma.simplices{isimp,1}.nodes(ind_nr(k0,1)), Gamma.simplices{isimp,1}.nodes(ind_nr(k0,2))];
        index_simp1(N1,1) = isimp;
    else
        N2 = N2+1;
        
        index_edges2(N2,1) = e;
        index_edges_corners2(N2,:) = [Gamma.simplices{isimp,1}.nodes(ind_nr(k0,1)), Gamma.simplices{isimp,1}.nodes(ind_nr(k0,2))];
        index_simp2(N2,1) = isimp;
    end
end
index_edges1 = index_edges1(1:N1,:);
index_edges_corners1 = index_edges_corners1(1:N1,:);
index_simp1 = index_simp1(1:N1,:);
index_edges2 = index_edges2(1:N2,:);
index_edges_corners2 = index_edges_corners2(1:N2,:);
index_simp2 = index_simp2(1:N2,:);


%% Get loop indices of corners 1
index1 = zeros(N1,1); 

index1(1,1) = index_edges_corners1(1,1); 
j=1;
for i=2:N1
    i0 = index_edges_corners1(j,2);
    % Find i0
    search = index_edges_corners1(:,1); 
    index_i0_in_search = find(search == i0); 
    
    j=index_i0_in_search;
    index1(i,1) = index_edges_corners1(j,1);
end
% Test 
i0 = index_edges_corners1(j,2); 
if(i0 == index1(1,1))
    fprintf('index1 okay, connected circle\n'); 
else
    fprintf('index1 false, non-circle\n'); 
end

    
% Get loop indices of corners 2
index2 = zeros(N2,1); 

index2(1,1) = index_edges_corners2(1,1); 
j=1;
for i=2:N2
    i0 = index_edges_corners2(j,2); 
    % Find i0
    search = index_edges_corners2(:,1); 
    index_i0_in_search = find(search == i0); 
    
    j=index_i0_in_search;
    index2(i,1) = index_edges_corners2(j,1); 
end
% Test 
i0 = index_edges_corners2(j,2); 
if(i0 == index2(1,1))
    fprintf('index2 okay, connected circle\n'); 
else
    fprintf('index2 false, non-circle\n'); 
end


%% Get cost fcn
Xc1 = zeros(N1,3); 
Xc2 = zeros(N2,3); 
for i=1:N1
    Xc1(i,:) = Gamma.X(index1(i,1),:); 
end
for i=1:N2
    Xc2(i,:) = Gamma.X(index2(i,1),:); 
end

C = zeros(N1,N2); 
for i=1:N1
    for j=1:N2
        C(i,j) = norm(Xc1(i,:)-Xc2(j,:)); 
    end
end

[Matching,Cost] = Hungarian(C);

%% Assignment and detection of free nodes
assign1 = [index1, zeros(N1,2)]; 
for i=1:N1
    search = Matching(i,:);
    k = find(search == 1);
    if(~isempty(k))
        assign1(i,2:3) = [k,index2(k)];
    end
end
assign2 = [index2, zeros(N2,2)]; 
for j=1:N2
    search = Matching(:,j); 
    k = find(search == 1);
    if(~isempty(k))
        assign2(j,2:3) = [k,index1(k)]; 
    end
end
    

%% Insert new points
insert1 = zeros(0,5); % isimp, edge, inode1,inode2,nr insert points
insert2 = zeros(0,5); 
Ninsert1 = 0; 
Ninsert2 = 0; 

for i=1:N1
    if(assign1(i,2) == 0)
        % Search in minus direction
        found = 0; 
        j=i-1;
        while(found == 0)
            if(j==0)
                j=N1;
            end
            if(assign1(j,2) > 0)
                found = 1;
            else
                j=j-1;
            end
        end
        Ninsert2 = Ninsert2+1; 
        insert2(Ninsert2,:) = [0,0,assign1(j,3),0,1];
        % Search in plus direction
        found = 0; 
        j=i+1; 
        while(found == 0)
            if(j==N1+1)
                j=1;
            end
            if(assign1(j,2) > 0)
                found = 1;
            else
                j=j+1;
            end
        end
        insert2(Ninsert2,4) = assign1(j,3); 
        % Find index of corresponding edge and simplex
        k = find(index_edges_corners2(:,1) == insert2(Ninsert2,3) & index_edges_corners2(:,2) == insert2(Ninsert2,4)); 
        if(isempty(k))
            k = find(index_edges_corners2(:,1) == insert2(Ninsert2,4) & index_edges_corners2(:,2) == insert2(Ninsert2,3)); 
        end
        insert2(Ninsert2,1:2) = [index_simp2(k,1), index_edges2(k,1)]; 
        
        % Look if simplex is already stored
        search = insert2(1:(Ninsert2-1),1); 
        k = find(search == insert2(Ninsert2,1)); 
        if(~isempty(k))
            insert2(Ninsert2,:) = []; 
            Ninsert2 = Ninsert2 - 1; 
            insert2(k,5) = insert2(k,5) + 1; 
        end
        
    end
end
for i=1:N2
    if(assign2(i,2) == 0)
        % Search in minus direction
        found = 0; 
        j=i-1;
        while(found == 0)
            if(j==0)
                j=N2;
            end
            if(assign2(j,2) > 0)
                found = 1;
            else
                j=j-1;
            end
        end
        Ninsert1 = Ninsert1+1; 
        insert1(Ninsert1,:) = [0,0,assign2(j,3),0,1];
        % Search in plus direction
        found = 0; 
        j=i+1; 
        while(found == 0)
            if(j==N2+1)
                j=1;
            end
            if(assign2(j,2) > 0)
                found = 1;
            else
                j=j+1;
            end
        end
        insert1(Ninsert1,4) = assign2(j,3); 
        % Find index of corresponding edge and simplex
        k = find(index_edges_corners1(:,1) == insert1(Ninsert1,3) & index_edges_corners1(:,2) == insert1(Ninsert1,4)); 
        if(isempty(k))
            k = find(index_edges_corners1(:,1) == insert1(Ninsert1,4) & index_edges_corners1(:,2) == insert1(Ninsert1,3)); 
        end
        insert1(Ninsert1,1:2) = [index_simp1(k,1), index_edges1(k,1)]; 
        
        % Look if simplex is already stored
        search = insert1(1:(Ninsert1-1),1); 
        k = find(search == insert1(Ninsert1,1)); 
        if(~isempty(k))
            insert1(Ninsert1,:) = []; 
            Ninsert1 = Ninsert1 - 1; 
            insert1(k,5) = insert1(k,5) + 1; 
        end
        
    end
end


[Gamma, index1, index_edges1, index_edges_corners1, index_simp1] = insert_new_points(Gamma,insert1,index1,index_edges1,index_simp1);
[Gamma, index2, index_edges2, index_edges_corners2, index_simp2] = insert_new_points(Gamma,insert2,index2,index_edges2,index_simp2);




%% Repeat Cost Function and Matching, final assignment
% Get cost fcn
N1 = size(index1,1); 
N2 = size(index2,1); 

Xc1 = zeros(N1,3); 
Xc2 = zeros(N2,3); 
for i=1:N1
    Xc1(i,:) = Gamma.X(index1(i,1),:); 
end
for i=1:N2
    Xc2(i,:) = Gamma.X(index2(i,1),:); 
end

C = zeros(N1,N2); 
for i=1:N1
    for j=1:N2
        C(i,j) = norm(Xc1(i,:)-Xc2(j,:)); 
    end
end

[Matching,Cost] = Hungarian(C);

if(Cost / N1 > 1)
    warning = 1;
else 
    warning = 0; 
end

if(warning == 0)
% Assignment and detection of free nodes
assign1 = [index1, zeros(N1,2)]; 
cost_vec = zeros(N1,1); 

for i=1:N1
    search = Matching(i,:);
    k = find(search == 1);
    if(~isempty(k))
        assign1(i,2:3) = [k,index2(k)];
        cost_vec(i,1) = C(i,k);
    end
end
assign2 = [index2, zeros(N2,2)]; 
for j=1:N2
    search = Matching(:,j); 
    k = find(search == 1);
    if(~isempty(k))
        assign2(j,2:3) = [k,index1(k)]; 
    end
end

%% Check if connections are valid
[value,istart] = min(cost_vec); 
jstart = assign1(istart,2); 

assign1_check = zeros(N1,3); 
assign2_check = zeros(N2,3); 
assign1_check(1,:) = [index1(istart), jstart, index2(jstart)]; 
assign2_check(1,:) = [index2(jstart), istart, index1(istart)]; 


search1 = index_edges_corners1(:,1); 
search2 = index_edges_corners2(:,2); 

i = find(search1 == index1(istart)); 
j = find(search2 == index2(jstart)); 

for count=2:N1
    inext = index_edges_corners1(i,2);
    i = find(search1 == inext);
    val = index_edges_corners1(i,1) ;
    
    assign1_check(count,1) = val;
    val_save = val; 
    
    jnext = index_edges_corners2(j,1);
    j = find(search2 == jnext);
    val = index_edges_corners2(j,2) ;
    
    
    assign1_check(count,2:3) = [find(index2 == jnext),val];
    
    assign2_check(count,:) = [val, find(index1 == inext), val_save]; 
    
end

failure = 0; 
for i=1:N1
    inode = assign1_check(i,1);
    search = assign1(:,1); 
    k = find(search == inode); 
    
    if(assign1(k,2) ~= assign1_check(i,2) || assign1(k,3) ~= assign1_check(i,3))
        failure = i;
    end
    inode = assign2_check(i,1); 
    search = assign2(:,1); 
    k = find(search == inode); 
    
    if(assign2(k,2) ~= assign2_check(i,2) || assign2(k,3) ~= assign2_check(i,3))
        failure = i;
    end
    
end



if(failure > 0)
    fprintf('Assignment orientation error! Use assign*_check instead of assign*\n'); 
    assign1 = assign1_check;
    assign2 = assign2_check; 
    pause
else
    fprintf('No failure in assignment orientation\n'); 
end




%% Get matching simplices
assign_simp = zeros(N1,6); % simp1,simp2,edge1,edge2,local index1,local index2

for i=1:N1
    isimp = index_simp1(i,1);
    iedge = index_edges1(i,1);
    
    i1 = index_edges_corners1(i,1);
    i2 = index_edges_corners1(i,2);
    
    % Search for edge/simplex of second group
    search = assign2(:,3);
    k1 = find(search == i1);
    j1 = assign2(k1,1);
    
    k2 = find(search == i2) ;
    j2 = assign2(k2,1);
    

    k = find(index_edges_corners2(:,1) == j1);
    Nfind = size(k,1);
    k0 = -1;
    for j=1:Nfind
        if(index_edges_corners2(k(j),2) == j2)
            k0 = k(j);
        end
    end
    if(Nfind==0 || k0==-1)
        k = find(index_edges_corners2(:,1) == j2);
        Nfind = size(k,1);
        for j=1:Nfind
            if(index_edges_corners2(k(j),2) == j1)
                k0 = k(j); 
            end
        end
    end
    
    

    assign_simp(i,:) = [isimp, index_simp2(k0,1), iedge, index_edges2(k0,1), 0, 0]; 
    
    % Get local index
    found = 0; 
    k=1;
    while(k<=3 && found==0)
        if(Gamma.simplices{isimp,1}.edges(k) == iedge)
            found = k; 
        end
        k=k+1;
    end
    assign_simp(i,5) = found; 
    
    found = 0; 
    k=1;
    while(k<=3 && found==0)
        if(Gamma.simplices{assign_simp(i,2),1}.edges(k) == assign_simp(i,4))
            found = k; 
        end
        k=k+1;
    end
    assign_simp(i,6) = found;     
end
else
    assign_simp = 0; 
end


Gamma_new = Gamma;
end





    
    