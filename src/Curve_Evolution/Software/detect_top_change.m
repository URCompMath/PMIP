function [split,merge,triple,boundary,info,mark] = detect_top_change(Gamma,Image,Surface,Prev,config,repeat)

%% Initializations / Preparations
X = Gamma.X; 
neigh = Gamma.neigh; 
index_info = Gamma.index_info; 
nr_curves = Gamma.nr_curves;
orient = Gamma.orient; 
I = Image.data; 
im_seg_method = config.method.flag; 


Nmain = size(Gamma.index_info,1);
Nsub  = size(Gamma.index_info,2);
mark  = zeros(Nmain,Nsub);



% Set size of neighborhood
size_nbh = 4; 

% Overlap method
if(im_seg_method == 2)
    overlap_flag = 1;
else
    overlap_flag = 0;
end

% Default output
split = zeros(0,4);    % node index 1, node index 2, main curve index, sub-curve index
merge = zeros(0,6);    % node index 1, main curve index 1 , sub-curve index 1, node index 2, main curve index 2 , sub-curve index 2, 
triple = zeros(0,6);   % node index 1, main curve index 1 , sub-curve index 1, node index 2, main curve index 2 , sub-curve index 2, 
boundary = zeros(0,4); % node index, main curve index, sub-curve index, boundary_index

nr_split = 0;
nr_merge = 0;
nr_triple = 0;
nr_boundary = 0;


%% Calculate grid size of underlying grid (later stored in cell info)
amax=config.top_check.amax;
amin=config.top_check.amin; 

norm_delta_X = sqrt(sum(Gamma.delta_X.*Gamma.delta_X,2)); 



a = 10*max(norm_delta_X); 
if(repeat>0 && config.dimension==2)  
    a = a * config.top_check.n_sub;
end
a = max([amin,a]); 
a = min([amax,a]); 
% fprintf('Top change: Grid size = %2.4f\n', a); 

 N0 = size(index_info,1); 

 %% Set info matrix / create sparse matrix info 
switch config.dimension
    case 2
        N = round(size(I,1)/a);
        M = round(size(I,2)/a);
        info = cell(N,M);
    case 3
        
 
        N = round((Surface.aminmax(2)-Surface.aminmax(1))/a) +1;
        M = round((Surface.aminmax(4)-Surface.aminmax(3))/a) +1;
        O = round((Surface.aminmax(6)-Surface.aminmax(5))/a) +1;
        info = cell(N,M,O);
end

%% Construct underlying grid to detect top changes efficiently
%% First loop: Set entries of info to [-1,-1,-1] if a grid point is in an array
% check also if a point touches the boundary (set entry of info to [-2,-2])
for ll=1:N0
    for l=1:nr_curves(ll,1)
        i = index_info(ll,l);
        count = 1; 
        
            

        while(i~=index_info(ll,l) || count == 1)
            switch config.dimension
                case 2
                    % Get index of grid array
                    k0 = min([max([1,ceil(X(i,2)/a)]),N]);
                    l0 = min([max([1,ceil(X(i,1)/a)]),M]);
                    
                    % Initialize info{k0,l0}
                    info{k0,l0}=[-1,-1,-1];
                    % Check if point is near a boundary
                    if(neigh(i,1)>0 && neigh(i,2)>0)  % no previous boundary point
                        boundary_index = is_boundary(X(i,:), size(I,2), size(I,1));
                        if(boundary_index > 0)
                            info{k0,l0}=[-2,-2,-2];  % mark field
                            
                            % Store new boundary point in boundary
                            nr_boundary = nr_boundary + 1;
                            boundary(nr_boundary,:) = [i,ll,l,boundary_index];
                        end
                    end
                case 3
                    % Get index of grid array
                    k0 = min([max([1,round((X(i,1)-Surface.aminmax(1))/a) +1;]),N]);
                    l0 = min([max([1,round((X(i,2)-Surface.aminmax(3))/a) +1;]),M]);
                    m0 = min([max([1,round((X(i,3)-Surface.aminmax(5))/a) +1;]),O]);
                    % Initialize info{k0,l0,m0}
                    info{k0,l0,m0}=[-1,-1,-1];
            end
            
                
                
            
            % Update i (go to next point along the current curve)
            i=neigh(i,2);
            if(i<0)
                i=index_info(ll,l);  %end loop
            end
            count = count + 1;
        end
    end
end

%% Second loop: Look for splitting/merging/etc. and store pairs
for ll=1:N0
    for l=1:nr_curves(ll,1)

        i = index_info(ll,l);
        count = 1; 
        while(i~=index_info(ll,l) || count == 1)

            switch config.dimension
                case 2
                    % Get index of grid array
                    k0 = min([max([1,ceil(X(i,2)/a)]),N]);
                    l0 = min([max([1,ceil(X(i,1)/a)]),M]);
                    info0 = info{k0,l0}(1);
                    Prev0 = in_Prev(X(i,:),Prev,config);
                    
                case 3
                    % Get index of grid array
                    k0 = min([max([1,round((X(i,1)-Surface.aminmax(1))/a) +1]),N]);
                    l0 = min([max([1,round((X(i,2)-Surface.aminmax(3))/a) +1]),M]);
                    m0 = min([max([1,round((X(i,3)-Surface.aminmax(5))/a) +1]),O]);
                    info0 = info{k0,l0,m0}(1);
                    Prev0 = in_Prev(X(i,:),Prev,config);

            end
                 

            % If grid array is not yet occupied store curve and node index
            % Do not perform top change if Prev(k0,l0)=1, i.e. if a top
            % change has happed there in the last timestep 
            if(info0<0 || Prev0>0 )
                switch config.dimension
                    case 2
                        info{k0,l0}=[ll,l,i];   % store curve and node index
                    case 3
                        info{k0,l0,m0}=[ll,l,i];
                end
            else 
                
                switch config.dimension
                    case 2
                        Help=info{k0,l0}; % get curve and node index stored in this array
                    case 3
                        Help=info{k0,l0,m0}; 
                end
                
                kk = Help(1);
                k = Help(2);
                j = Help(3);

                
                if(kk==ll && k==l)  % nodes belong to the same (main and) sub-curve --> splitting of the sub-curve 
                    % perform splitting only if non-boundary nodes are involved
                    if(neigh(i,1)>0 && neigh(i,2)>0 && neigh(j,1)>0 && neigh(j,2)>0) 
                        if(neighbor(i,j,neigh,size_nbh)==0)  %no neighbor points

                            [i1,j1] = local_check(X,i,j,neigh);  
                            if(i1>0 && j1>0)
                                [i1,j1] = sort_indices(i1,j1,ll,l,index_info,neigh);
                                nr_split = nr_split + 1;
                                split(nr_split,:) = [i1,j1,ll,l];
                                
                                % jump 2 points forward (or go to end of curve)
                                i_latest = find_latest(neigh,i,i1,index_info(ll,l));  % fcn find_latest modified 2013-07-20
                                
                                if(neigh(i_latest,2)~=index_info(ll,l) && neigh(i_latest,2)>0)
                                    if(neigh(neigh(i_latest,2),2)~= index_info(ll,l) && neigh(neigh(i_latest,2),2)>0)
                                        i=neigh(neigh(i_latest,2),2);
                                    else
                                        i=neigh(index_info(ll,l),1);
                                        % If index_info is a neighbor and
                                        % already in the list, delete entry
                                        if(i_in_list(boundary,split,merge,triple,index_info(ll,l)))
                                            split(nr_split,:)=[];
                                            nr_split = nr_split - 1;
                                        end
                                   end
                                else  % go to end of curve
                                    i=neigh(index_info(ll,l),1);

                                    % If index_info is a neighbor and
                                    % already in the list, delete entry
                                    if(i_in_list(boundary,split,merge,triple,index_info(ll,l)))
                                        split(nr_split,:)=[];
                                        nr_split = nr_split - 1;
                                    end
                                end
                                
                            end
                        end
                    end
                    
                else  % nodes do not belong to the same sub-curve --> merging of the two sub-curves or creation of triple junctions
                    % perform top change only if non-boundary nodes are involved

                    size_bd = 3;  

                       if(near_boundary(i,neigh,size_bd)==0 && near_boundary(j,neigh,size_bd)==0)
                           [i1,j1] = local_check(X,i,j,neigh);
                           if(i1>0 && j1>0)
                               if((orient(kk,1)==orient(ll,1) && orient(kk,2)==orient(ll,2)) || (orient(kk,1)==orient(ll,2) && orient(kk,2)==orient(ll,1)))
                                   % Merge curves
                                   nr_merge = nr_merge + 1;
                                   merge(nr_merge,:) = [i1,ll,l,j1,kk,k];
                                   flag = 1;
                               else
                                   if(~overlap_flag)
                                       % Triple junction
                                       nr_triple = nr_triple + 1;
                                       triple(nr_triple,:) = [i1,ll,l,j1,kk,k];

                                   end
                                   flag = 2;
                               end
                               % jump 2 points forward (or go to end of curve)
                               i_latest = find_latest(neigh,i,i1,index_info(ll,l));  % fcn find_latest modified 2013-07-20
                               if(neigh(i_latest,2)~=index_info(ll,l) && neigh(i_latest,2)>0)
                                   if(neigh(neigh(i_latest,2),2)~= index_info(ll,l) && neigh(neigh(i_latest,2),2)>0)
                                       i=neigh(neigh(i_latest,2),2);
                                   else
                                       i=neigh(index_info(ll,l),1);
                                       % If index_info is a neighbor and
                                       % already in the list, delete entry
                                       if(i_in_list(boundary,split,merge,triple,index_info(ll,l)))
                                           if(flag == 1)
                                               merge(nr_merge,:)=[];
                                               nr_merge = nr_merge - 1;
                                           else
                                               triple(nr_triple,:)=[];
                                               nr_triple = nr_triple - 1;
                                           end
                                       end
                                       
                                   end
                               else  % go to end of curve
                                   i=neigh(index_info(ll,l),1);
                                   % If index_info is a neighbor and
                                   % already in the list, delete entry
                                   if(i_in_list(boundary,split,merge,triple,index_info(ll,l)))
                                       if(flag == 1)
                                           merge(nr_merge,:)=[];
                                           nr_merge = nr_merge - 1;
                                       else
                                           triple(nr_triple,:)=[];
                                           nr_triple = nr_triple - 1;
                                       end
                                   end
                               end
                           end
                       end
                       

                end
            end

            if(i>0)
               i=neigh(i,2);
            end
            count = count + 1;
            if(i<0)
                i=index_info(ll,l);
            end
            
        end
    end
end


%% Check if a node appears two times
n_boundary = size(boundary,1); 
n_split = size(split,1); 
n_merge = size(merge,1); 
n_triple = size(triple,1); 

% loop over boundary indices
for ib=1:n_boundary
    i=boundary(ib,1);
    ll=boundary(ib,2); 
    l=boundary(ib,3); 
    
    if(i>0)
        collect = zeros(0,1);
        count = 0; 
        for j=(ib+1):n_boundary
            if(i==boundary(j,1) && ll==boundary(j,2) && l==boundary(j,3))
                count = count + 1; 
                collect(count,1)=j;
            end
        end
        for j=1:count
            boundary(collect(j,1),:)=zeros(1,4); 
        end
        collect = zeros(0,1);
        count = 0; 
        for j=1:n_split
            if((i==split(j,1) || i==split(j,2)) && ll==split(j,3) && l==split(j,4))
                count = count + 1; 
                collect(count,1)=j;
            end
        end
        for j=1:count
            split(collect(j,1),:)=zeros(1,4); 
        end
        collect = zeros(0,1);
        count = 0; 
        for j=1:n_merge
            if((i==merge(j,1) && ll==merge(j,2) && l==merge(j,3)) || (i==merge(j,4) && ll==merge(j,5) && l==merge(j,6)))
                count = count + 1; 
                collect(count,1)=j;
            end
        end
        for j=1:count           
            merge(collect(j,1),:)=zeros(1,6); 
        end
        collect = zeros(0,1);
        count = 0; 
        for j=1:n_triple
            if((i==triple(j,1) && ll==triple(j,2) && l==triple(j,3)) || (i==triple(j,4) && ll==triple(j,5) && l==triple(j,6)))
                count = count + 1; 
                collect(count,1)=j;
            end
        end
        for j=1:count           
            triple(collect(j,1),:)=zeros(1,6); 
        end
    end    
end

% loop over split indices
for is=1:n_split
    i = split(is,1); 
    j1 = split(is,2);
    ll = split(is,3); 
    l = split(is,4); 
    
    if(i>0)
        collect = zeros(0,1);
        count = 0; 
        for j=(is+1):n_split
            if(((i==split(j,1) && j1==split(j,2)) || (i==split(j,2) && j1==split(j,1))) && ll==split(j,3) && l==split(j,4))
                count = count + 1; 
                collect(count,1)=j;
            end
        end
        for j=1:count
            split(collect(j,1),:)=zeros(1,4); 
        end
        collect = zeros(0,1);
        count = 0; 
        for j=1:n_merge
            if((i==merge(j,1) && ll==merge(j,2) && l==merge(j,3)) || (i==merge(j,4) && ll==merge(j,5) && l==merge(j,6)) || (j1==merge(j,1) && ll==merge(j,2) && l==merge(j,3)) || (j1==merge(j,4) && ll==merge(j,5) && l==merge(j,6)) )
                count = count + 1; 
                collect(count,1)=j;
            end
        end
        for j=1:count           
            merge(collect(j,1),:)=zeros(1,6); 
        end
        collect = zeros(0,1);
        count = 0; 
        for j=1:n_triple
            if((i==triple(j,1) && ll==triple(j,2) && l==triple(j,3)) || (i==triple(j,4) && ll==triple(j,5) && l==triple(j,6)) || (j1==triple(j,1) && ll==triple(j,2) && l==triple(j,3)) || (j1==triple(j,4) && ll==triple(j,5) && l==triple(j,6) ))
                count = count + 1; 
                collect(count,1)=j;
            end
        end
        for j=1:count           
            triple(collect(j,1),:)=zeros(1,6); 
        end
    end    
end
% loop over merge indices
for im=1:n_merge
    i = merge(im,1); 
    ll = merge(im,2); 
    l = merge(im,3); 
    j1 = merge(im,4); 
    kk = merge(im,5); 
    k = merge(im,6); 
    
    if(i>0)
        collect = zeros(0,1);
        count = 0; 
        for j=(im+1):n_merge
            if((i==merge(j,1) && ll==merge(j,2) && l==merge(j,3)) || (i==merge(j,4) && ll==merge(j,5) && l==merge(j,6)) || (j1==merge(j,1) && kk==merge(j,2) && k==merge(j,3)) || (j1==merge(j,4) && kk==merge(j,5) && k==merge(j,6)) )
                count = count + 1; 
                collect(count,1)=j;
            end
        end
        for j=1:count           
            merge(collect(j,1),:)=zeros(1,6); 
        end
        collect = zeros(0,1);
        count = 0; 
        for j=1:n_triple
            if((i==triple(j,1) && ll==triple(j,2) && l==triple(j,3)) || (i==triple(j,4) && ll==triple(j,5) && l==triple(j,6)) || (j1==triple(j,1) && kk==triple(j,2) && k==triple(j,3)) || (j1==triple(j,4) && kk==triple(j,5) && k==triple(j,6) ))
                count = count + 1; 
                collect(count,1)=j;
            end
        end
        for j=1:count           
            triple(collect(j,1),:)=zeros(1,6); 
        end
    end    
end
% loop over triple indices
for it=1:n_triple
    i = triple(it,1); 
    ll = triple(it,2); 
    l = triple(it,3); 
    j1 = triple(it,4); 
    kk = triple(it,5); 
    k = triple(it,6); 
    
    if(i>0)
        collect = zeros(0,1);
        count = 0; 
        for j=(it+1):n_triple
            if((i==triple(j,1) && ll==triple(j,2) && l==triple(j,3)) || (i==triple(j,4) && ll==triple(j,5) && l==triple(j,6)) || (j1==triple(j,1) && kk==triple(j,2) && k==triple(j,3)) || (j1==triple(j,4) && kk==triple(j,5) && k==triple(j,6) ))
                count = count + 1; 
                collect(count,1)=j;
            end
        end
        for j=1:count           
            triple(collect(j,1),:)=zeros(1,6); 
        end
    end    
end

%% Shift entries
collect = zeros(0,1); 
count = 0; 
for i=1:n_boundary
    if(boundary(i,1)==0)
        count = count +1; 
        collect(count,1)=i;
    end
end
for i=1:count
    l = collect(i,1)-i+1; 
    boundary(l:(n_boundary-1),:)=boundary((l+1):n_boundary,:);
    boundary(n_boundary,:) = zeros(1,4); 
end
n_new = n_boundary - count; 
boundary((n_new+1):n_boundary,:)=[]; 

collect = zeros(0,1); 
count = 0; 
for i=1:n_split
    if(split(i,1)==0)
        count = count +1; 
        collect(count,1)=i;
    end
end
for i=1:count
    l = collect(i,1)-i+1; 
    split(l:(n_split-1),:)=split((l+1):n_split,:);
    split(n_split,:) = zeros(1,4); 
end
n_new = n_split - count; 
split((n_new+1):n_split,:)=[]; 

collect = zeros(0,1); 
count = 0; 
for i=1:n_merge
    if(merge(i,1)==0)
        count = count +1; 
        collect(count,1)=i;
    end
end
for i=1:count
    l = collect(i,1)-i+1; 
    merge(l:(n_merge-1),:)=merge((l+1):n_merge,:);
    merge(n_merge,:) = zeros(1,6); 
end
n_new = n_merge - count; 
merge((n_new+1):n_merge,:)=[]; 

collect = zeros(0,1); 
count = 0; 
for i=1:n_triple
    if(triple(i,1)==0)
        count = count +1; 
        collect(count,1)=i;
    end
end
for i=1:count
    l = collect(i,1)-i+1; 
    triple(l:(n_triple-1),:)=triple((l+1):n_triple,:);
    triple(n_triple,:) = zeros(1,6); 
end
n_new = n_triple - count; 
triple((n_new+1):n_triple,:)=[]; 


% delete double entries  
[split, merge, triple] = delete_double(split,merge,triple);

end


function out = neighbor(i,j,neigh,size_nbh)
if((neigh(i,1)<0 && neigh(j,2)<0) || (neigh(i,2)<0 && neigh(j,1)<0))
    out = 1;  % prevent splitting in curves with one curve has length 0 
else
    
    size_nbh=round(size_nbh);

    info1 = zeros(size_nbh,1); 
    info2 = zeros(size_nbh,1); 

    stop = 0; 
    k=0; 
    kk = i;
    while(k<= size_nbh && stop==0)
        kk = neigh(kk,1); 
        if(kk>0)
            k = k+1;
            info1(k,1)=kk;
        else
            stop=1;
        end
    end
    info1=info1(1:k,1);

    stop = 0;
    k=0;
    kk=i;
    while(k<= size_nbh && stop==0)
        kk = neigh(kk,2); 
        if(kk>0)
            k = k+1;
            info2(k,1)=kk;
        else
            stop=1;
        end
    end
    info2=info2(1:k,1);
    info = [info1;info2];



    % Default output
    out = 0;

    % Set out to 1, if j is in info
    k=1;
    while(k<=size(info,1))
        if(j==info(k))
            out = 1;
            k = size(info,1);
        end
        k = k+1;
    end
end        

end


function [i0,j0] = local_check(X,i,j,neigh)

I = zeros(1,5); 
J = zeros(1,5);


% index set I
K=0; 
if(neigh(i,1)>0 && neigh(i,2)>0 && neigh(neigh(i,1),1)>0 && neigh(neigh(i,2),2)>0)
    I(1)=i;
    K=1;
end

% loop in + direction
stop=0;
k=1;
ii=i;
while(stop==0 && k<=5)
    ii = neigh(ii,2);
    if(ii>0 && neigh(ii,1)>0 && neigh(ii,2)>0 && neigh(neigh(ii,1),1)>0 && neigh(neigh(ii,2),2)>0)
        I(K+k)=ii;
        k=k+1;
    else
        if(ii<0 || neigh(ii,2)<0)
            stop=1;
        else
            if(neigh(neigh(ii,2),2)<0)
                stop=1;
            end
        end
    end
end
K=K+k-1;

% loop in - direction
stop=0; 
k=1; 
ii = i; 
while(stop==0 && k<=5)
    ii = neigh(ii,1);
    if(ii>0 && neigh(ii,1)>0 && neigh(ii,2)>0 && neigh(neigh(ii,1),1)>0 && neigh(neigh(ii,2),2)>0 )
        I(K+k)=ii;
        k=k+1;
    else
        if(ii<0 || neigh(ii,1)<0)
        stop=1;
        else
            if(neigh(neigh(ii,1),1)<0)
                stop=1;
            end
        end
    end
end
K=K+k-1;


I=I(1,1:K); 

% index set J
L=0; 
if(neigh(j,1)>0 && neigh(j,2)>0 && neigh(neigh(j,1),1)>0 && neigh(neigh(j,2),2)>0 )
    J(1)=j;
    L=1;
end

% loop in + direction
stop=0;
k=1;
ii=j;
while(stop==0 && k<=5)
    ii = neigh(ii,2);
    if(ii>0 && neigh(ii,1)>0 && neigh(ii,2)>0 && neigh(neigh(ii,1),1)>0 && neigh(neigh(ii,2),2)>0 )
        J(L+k)=ii;
        k=k+1;
    else
        if(ii<0 || neigh(ii,2)<0)
            stop=1;
        else
            if(neigh(neigh(ii,2),2)<0)
                stop=1;
            end
        end
    end
end
L=L+k-1;

% loop in - direction
stop=0; 
k=1; 
ii = j; 
while(stop==0 && k<=5)
    ii = neigh(ii,1);
    if(ii>0 && neigh(ii,1)>0 && neigh(ii,2)>0 && neigh(neigh(ii,1),1)>0 && neigh(neigh(ii,2),2)>0)
        J(L+k)=ii;
        k=k+1;
    else
        if(ii<0 || neigh(ii,1)<0)
            stop=1;
        else
            if(neigh(neigh(ii,1),1)<0)
                stop=1;
            end
        end
    end
end
L=L+k-1;


J=J(1,1:L); 


% Compute pair with smallest distance
% Initialize dist, i0, j0:
dist = Inf; 
i0 = -1;
j0 = -1;

% Loop to find pair (i0,j0) of smallest distance
for k=1:K
    for l=1:L
        dist0 = norm(X(I(k),:)-X(J(l),:));
        if(dist0<dist && I(k)~=J(l))
            dist = dist0;
            i0 = I(k);
            j0 = J(l);
        end
    end
end


end


function [i_out,j_out] = sort_indices(i_in,j_in,gamma_nr,curve_nr,index_info,neigh)

count_i = 0;
count_j = 0;
i= index_info(gamma_nr,curve_nr);
j= index_info(gamma_nr,curve_nr);

while(i ~= i_in)
    i=neigh(i,2);
    count_i = count_i + 1;
end

while(j ~= j_in)
    j = neigh(j,2);
    count_j = count_j + 1;
end

if(count_i >= count_j)
    i_out = i_in;
    j_out = j_in;
else 
    i_out = j_in;
    j_out = i_in;
end


    
end



function out = near_boundary(i,neigh,size_bd)
out  = 0;
stop = 0;
k    = 1;
i0=i;

while(k<=size_bd && stop==0)
    ii=neigh(i,1);
    if(ii>0)
        i=ii;
        k=k+1;
    else
        if(ii==-1)
            stop=-1;
        else
            stop = 1;
        end
    end
end
if(stop==1)
    out = 1;
else
    k=1; 
    i=i0;
    while(k<=size_bd && stop==0)
        ii=neigh(i,2);
        if(ii>0)
            i=ii;
            k=k+1;
        else
            if(ii==-1)
                stop=-1;
            else
                stop = 1;
            end
        end
    end
    if(stop==1)
        out = 1;
    end
end



end

function out = i_in_list(boundary,split,merge,triple,index)
out = 0; 
n_boundary = size(boundary,1); 
n_split = size(split,1);
n_merge = size(merge,1); 
n_triple = size(triple,1); 

i = 1; 
while(i <= n_boundary && out==0) 
    if(boundary(i,1)==index)
        out = 1;
    end
    i=i+1;
end
i=1; 
while(i <= n_split && out==0) 
    if(split(i,1)==index || split(i,2)==index)
        out = 1;
    end
    i=i+1;
end
i=1; 
while(i <= n_merge && out==0) 
    if(merge(i,1)==index || merge(i,4)==index)
        out = 1;
    end
    i=i+1;
end
i=1; 
while(i <= n_triple && out==0) 
    if(triple(i,1)==index || triple(i,4)==index)
        out = 1;
    end
    i=i+1;
end
end


function out = in_Prev(p,Prev,config)
n = size(Prev,1);


found = 0; 
j=1; 
while(j<=n && found==0)
    if(norm(Prev(j,:)-p)<config.top_check.a_Prev)
        found = 1; 
    end
    j=j+1;
end

out = found;
end


function out = find_latest(neigh,i0,i1,istart)
out = i0;

i=i0;
j=1;
while(i ~= i0 || j==1)
    i = neigh(i,2) ;
    j = j+1; 
    
    % stop loop
    if(i<0 || i==istart)  
        i=i0;
    else
        if(i==i1)
            out = i1;
            i=i0; 
        end
    end
end
end



function [split0, merge0, triple0] = delete_double(split,merge,triple)

% split
% split: Nx4    % node index 1, node index 2, main curve index, sub-curve index
N = size(split,1); 
delete = zeros(N,1); 
for i=1:N
    if(delete(i)==0)
        main1 = split(i,3); 
        sub1  = split(i,4); 
        
        for j=(i+1):N
            if(delete(j)==0)
                if(split(j,3)==main1 && split(j,4)==sub1)
                    delete(j)=1; 
                end
            end
        end
    end
end

split0 = split;
k = find(delete==1); 
split0(k,:) = []; 

% merge
% merge: Nx6    % node index 1, main curve index 1 , sub-curve index 1, node index 2, main curve index 2 , sub-curve index 2, 
N = size(merge,1); 
delete = zeros(N,1); 
for i=1:N
    if(delete(i)==0)
        main1 = merge(i,2); 
        sub1  = merge(i,3); 
        main2 = merge(i,5); 
        sub2  = merge(i,6); 
        
        for j=(i+1):N
            if(delete(j)==0)
                if(merge(j,2)==main1 && merge(j,3)==sub1 && merge(j,5)==main2 && merge(j,6)==sub2)
                    delete(j)=1; 
                end
            end
        end
    end
end 

merge0 = merge;
k = find(delete==1); 
merge0(k,:)=[]; 

% triple
% triple: Nx6  % node index 1, main curve index 1 , sub-curve index 1, node index 2, main curve index 2 , sub-curve index 2, 
N = size(triple,1); 
delete = zeros(N,1); 
for i=1:N
    if(delete(i)==0)
        main1 = triple(i,2); 
        sub1  = triple(i,3); 
        main2 = triple(i,5); 
        sub2  = triple(i,6); 
        
        for j=(i+1):N
            if(delete(j)==0)
                if(triple(j,2)==main1 && triple(j,3)==sub1 && triple(j,5)==main2 && triple(j,6)==sub2)
                    delete(j)=1; 
                end
            end
        end
    end
end 

triple0 = triple;
k = find(delete==1); 
triple0(k,:)=[]; 
end
