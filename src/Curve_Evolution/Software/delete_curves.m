function Gamma_new=delete_curves(Gamma)
%% Preparations
% Preparations
X            = Gamma.X;
neigh        = Gamma.neigh;
length_info  = Gamma.length;
index_info   = Gamma.index_info;
nr_curves    = Gamma.nr_curves;
orient       = Gamma.orient;
Lambda       = Gamma.Lambda;
mark         = Gamma.mark;
age          = Gamma.age; 

dim = size(X,2); 
if(dim==3)
    closest_simp = Gamma.closest_simp;
end


%% Delete curves
N0 = size(index_info,1);

mark_list = zeros(0,2);
mark_list_count = 0;

for ll=1:N0
    for l=1:nr_curves(ll,1)
        if(mark(ll,l)==1)
            mark_list_count = mark_list_count + 1;
            mark_list(mark_list_count,:)=[ll,l];
        end
    end
end


%% Delete curves one after another
for mark_count = 1:mark_list_count    
    ll = mark_list(mark_count,1);
    l  = mark_list(mark_count,2);
%     nr_points_deleted =0;
    
    if(ll>0)  % mark_list row can be set to -1 --> curve is previously deleted
        count = 0;
        J  = size(X,1);
        info = zeros(J,1);
        
        % Additional curve indices to be deleted
        add_delete = zeros(2,size(index_info,1)*size(index_info,2));
        add_delete_count = 0;
        row_delete_Lambda = zeros(1,1);
        row_delete_Lambda_count = 0;
        
        % New curves (new triple junctions)
        new_Lambda_info = zeros(0,4);
        new_Lambda_info_count = 0;
        
        i=index_info(ll,l);
        j = 1;
        istriple = 0;
        if(neigh(i,1)==-2)
            istriple=1;
            i_search = i;
        end
        if(neigh(i,1)<0)
            i_start = i;
        end
        while(i~=index_info(ll,l) || j == 1)  % get i_end
            i_end = i;
            j = j + 1;
            i=neigh(i,2);
            if(i<0)
                i=index_info(ll,l);
            end
        end
        
        
        if(neigh(i_end,2)==-3)
            istriple=1;
            if(neigh(i_start,1)>-2)  % boundary contact at i_start (-1), triple junction at i_end (-3)
                i_search = i_end;
            end
        end
        
        if(istriple)
            [row1,col1]=find(Lambda(:,1:3)==i_search);
            points_start = Lambda(row1,:);
            points_start(:,[col1,col1+3])=[];
            if(neigh(i_start,1)==-2 && neigh(i_end,2)==-3 )  % two triple points --> close curves, i_search = i_start
                [row2,col2]=find(Lambda(:,1:3) == i_end);
                points_end = Lambda(row2,:);
                points_end(:,[col2,col2+3])=[];
                
                % Now: in points_start --> other indices of nodes at
                % triple point belonging to i_start
                % in points_end --> other indices of nodes at triple
                % point belonging to i_end, indices in points_start
                % need not be connected to indices in points_end (e.g.
                % if a network of triple points exists)
                
                
                found_match = zeros(2,1);
                for j=1:2
                    
                    % Check if entry in index_info needs to be deleted
                    i0 = points_start(j);
                    main_curve_ind_i0 = points_start(j+2);
                   
                    
                    
                    for jj=3:4
                        if(points_end(jj)==main_curve_ind_i0)
                            found_match(j) = 1;
                            iN=points_end(jj-2);
                        end
                    end
                    
                    if(found_match(j))
                        if(neigh(i0,1)==-2)
                            i_run = i0;
                            iter = 1;
                            iN2 = i0;
                            while(i_run ~= i0  || iter == 1)
                                iN2 = i_run;
                                i_run = neigh(i_run,2); % run forwards
                                if(i_run<0)
                                    i_run = i0;
                                end
                                iter = iter + 1;
                                
                            end
                            if(iN2~=iN)  % delete index_info entry with i0, 2 subcurves --> close open hole at i0,iN --> one subcurve
                                [k0,l0]=find(index_info==i0);
                                add_delete_count = add_delete_count  +1;
                                add_delete(:,add_delete_count)=[k0;l0];
                            end
                            
                            % Change neighbor information, open curve (with triple
                            % junction(s) --> close curve at i0,iN)
                            neigh(i0,1)=iN;
                            neigh(iN,2)=i0;
                            
                        else % if neigh(i0,2)==-3
                            
                            i_run = i0;
                            iter = 1;
                            iN2 = i0;
                            while(i_run ~= i0 || iter == 1)
                                iN2 = i_run;
                                i_run = neigh(i_run,1); % run backwards
                                if(i_run<0)
                                    i_run = i0;
                                end
                                iter = iter + 1;
                            end

                            if(iN2~=iN) % delete index_info entry with iN, 2 subcurves --> close open hole at i0,iN --> one subcurve
                                [k0,l0]=find(index_info==iN);
                                add_delete_count = add_delete_count + 1;
                                add_delete(:,add_delete_count)=[k0;l0];
                            end
                            
                            % Change neighbor information,  open curve (with triple
                            % junction(s) --> close curve at i0,iN)
                            neigh(i0,2)=iN;
                            neigh(iN,1)=i0;
                            
                        end
                        
                    end
                    
                    
                end
                
                if(found_match(1)+found_match(2)==0) % 4 different curves meet at a quadrupel junction after deletion of marked curve
%                   fprintf('found_match(1)+found_match(2)==0, quadrupel junction with 4 diff curves\n'); 

                    % as the 4 curves are different (no matching of the
                    % same main curve) the 4 phases are different (can be
                    % proved by contradiction, see also quadrupel.pdf)
                    
                    % Curve points_start(3) seperates phase 1 and phase 2
                    % Curve points_end(3)   seperates phase 3 and phase 4
                    
                    
                    
                    new_Lambda_info_count = new_Lambda_info_count + 2;
                    phase1 = orient(points_start(3),1);
                    phase2 = orient(points_start(3),2);
                    phase3 = orient(points_end(3),1);
                    phase4 = orient(points_end(3),2);
                    
                    % Checking for phase 1 == phase 3 (example), is used to
                    % find out, if curve points_start(3) and curve
                    % points_end(3) lie next to each other or opposite
                    
                    if(phase1==phase3 || phase1==phase4 || phase2==phase3 || phase2==phase4)
                        % Curve start(3) and end(3) lie next to each other
                        % new_Lambda_info =[node1,node2,kmain1,kmain2]
                        new_Lambda_info(new_Lambda_info_count-1,:) = [points_start(1), points_end(1), points_start(3), points_end(3)];
                        new_Lambda_info(new_Lambda_info_count,:)   = [points_start(2), points_end(2), points_start(4), points_end(4)];
                    else
                        % Curve start(3) and end(4) lie next to each other
                        new_Lambda_info(new_Lambda_info_count-1,:) = [points_start(1), points_end(2), points_start(3), points_end(4)];
                        new_Lambda_info(new_Lambda_info_count,:)   = [points_start(2), points_end(1), points_start(4), points_end(3)];
                    end
                    
                end
                % Now: new triple junctions are stored in new_Lambda_info
                
                % Delete rows in Lambda corresponding to the former 2
                % triple junctions --> store information in
                % row_delete_Lambda
                row_delete_Lambda_count = row_delete_Lambda_count + 1;
                row_delete_Lambda(row_delete_Lambda_count,1)=row1;
                row_delete_Lambda_count = row_delete_Lambda_count + 1;
                row_delete_Lambda(row_delete_Lambda_count,1)=row2;
                
            
            else % neigh(i_end,2)==-1 || neigh(i_start,1)==-1 --> Curve connects one triple point and one boundary point
                row_delete_Lambda_count = row_delete_Lambda_count + 1;
                row_delete_Lambda(row_delete_Lambda_count,1)=row1;
                
                
                if(i_search == i_start)  % i_end boundary point, i_start triple point
                    X0 = X(i_end,:);
                    
                else  % i_start boundary point, i_end triple point (i_search == i_end)
                    X0 = X(i_start,:);
                    
                end
                % Project to boundary and change neighbor
                % information to -1
                for j=1:2
                    X(points_start(j),:)=X0;
                    if(neigh(points_start(j),1)==-2)
                        neigh(points_start(j),1)=-1;
                    else
                        if(neigh(points_start(j),2)==-3)
                            neigh(points_start(j),2)=-1;
                        end
                    end
                end
                
            end
            
        end
        i=index_info(ll,l);
        j=1;
        
        %% Curve (ll,l) has to be deleted, collect all indices of nodes
        % belonging to the curve and store the indices in info
        
        while(i~=index_info(ll,l) || j == 1)  % collect all indices i to in info to be deleted
            % Set info(count) to i
            info(count+1)=i;
            count = count +1;
            i_end = i;
            j = j + 1;
            i=neigh(i,2);
            if(i<0)
                i=index_info(ll,l);
            end
        end
        
        nr_points_deleted = count;
        
        
        %% Create new triple junction and new curve if new_Lambda_info_count > 0
        for count = 1:(new_Lambda_info_count/2)
            % Get nodes and curve info
            tp1 = new_Lambda_info(2*count-1,:);
            tp2 = new_Lambda_info(2*count,:);
            
            i11 = tp1(1);
            i12 = tp1(2);
            k11 = tp1(3);
            k12 = tp1(4);
            i21 = tp2(1);
            i22 = tp2(2);
            k21 = tp2(3);
            k22 = tp2(4);
            
            
            % Project nodes, normally: center between i11,i12 and i21,i22 resp.
            % However centers will be the same as X(i11,:)=X(i21,:) and
            % X(i12,:)=X(i22,:) before the deletion of the curve, so also include
            % neighbor points, to avoid the new curve to have length 0!
            
            if(neigh(i11,1)==-2)
                n11 = neigh(i11,2);
            else
                n11 = neigh(i11,1);
            end
            if(neigh(i12,1)==-2)
                n12 = neigh(i12,2);
            else
                n12 = neigh(i12,1);
            end
            if(neigh(i21,1)==-2)
                n21 = neigh(i21,2);
            else
                n21 = neigh(i21,1);
            end
            if(neigh(i22,1)==-2)
                n22 = neigh(i22,2);
            else
                n22 = neigh(i22,1);
            end
            
            
            X01 = (X(i11,:)+X(n11,:)+X(i12,:)+X(n12,:))/4;
            X02 = (X(i21,:)+X(n21,:)+X(i22,:)+X(n22,:))/4;
            
            X(i11,:) = X01;
            X(i12,:) = X01;
            X(i21,:) = X02;
            X(i22,:) = X02;
            
            % Project to surface and update closest simp if dim==3
            if(dim==3)
                [closest_simp(i11,1), X(i11,:)] = get_closest_simp_use_parents(X(i11,:),closest_simp(i11,1),Surface); 
                [closest_simp(i12,1), X(i12,:)] = get_closest_simp_use_parents(X(i12,:),closest_simp(i12,1),Surface); 
                [closest_simp(i21,1), X(i21,:)] = get_closest_simp_use_parents(X(i21,:),closest_simp(i21,1),Surface); 
                [closest_simp(i22,1), X(i22,:)] = get_closest_simp_use_parents(X(i22,:),closest_simp(i22,1),Surface); 
            end
            
            
            % Create 3 new points
            J=size(X,1);
            J=J+3;
            X(J-2,:) = X01;
            X(J-1,:) = (X01+X02)/2;
            X(J,:)   = X02;
            
            % Project to surface and update closest simp if dim==3
            if(dim==3)
                [closest_simp(J-2,1),X(J-2,:)] = get_closest_simp_use_parents(X(J-2,:),closest_simp(i11,1),Surface); 
                [closest_simp(J-1,1),X(J-1,:)] = get_closest_simp_use_parents(X(J-1,:),closest_simp(J-2,1),Surface); 
                [closest_simp(J  ,1),X(J  ,:)] = get_closest_simp_use_parents(X(J  ,:),closest_simp(i21,1),Surface); 
            end    
            
            % Get orientation
            phase11a = orient(k11,1);
            phase11b = orient(k11,2);
            phase12a = orient(k12,1);
            phase12b = orient(k12,2);
            
            if(phase11a == phase12a)
                common_phase = phase11a;
                phase_new1   = phase11b;
                phase_new2   = phase12b;
            else
                if(phase11a == phase12b)
                    common_phase = phase11a;
                    phase_new1   = phase11b;
                    phase_new2   = phase12a;
                else
                    if(phase11b == phase12a)
                        common_phase = phase11b;
                        phase_new1   = phase11a;
                        phase_new2   = phase12b;
                    else
                        common_phase = phase11b;
                        phase_new1   = phase11a;
                        phase_new2   = phase12a;
                    end
                end
            end
            
            
            % Default new orient line
            new_orient_line = [phase_new1, phase_new2];
            if((neigh(i11,1)==-2 && common_phase == phase11b) || (neigh(i11,2)==-3 && common_phase == phase11a))
                new_orient_line = [phase_new2, phase_new1];
            end
            
            % Update index_info etc. if orient line already exists
            
            N0 = size(orient,1);
            change = 0;
            found = 0;
            for k=1:N0
                if(orient(k,1)==new_orient_line(1) && orient(k,2) == new_orient_line(2))
                    found = 1;
                    kmain = k;
                    nr_curves(k,1)=nr_curves(k,1)+1;
                    index_info(k,nr_curves(k,1))=J-2;
                    length_info(k,nr_curves(k,1))=norm(X(J,:)-X(J-2,:));
                    age(k,nr_curves(k,1))=1; 
                else
                    if(orient(k,1)==new_orient_line(2) && orient(k,2) == new_orient_line(1))
                        found = 1;
                        kmain = k;
                        nr_curves(k,1)=nr_curves(k,1)+1;
                        index_info(k,nr_curves(k,1))=J; 
                        length_info(k,nr_curves(k,1))=norm(X(J,:)-X(J-2,:));
                        age(k,nr_curves(k,1))=1; 
                        change = 1;
                    end
                end
            end
            
            % Create new main curve if orient line does not alreay exist
            if(found == 0)
                kmain = N0+1;
                nr_curves(kmain,1)=1;
                index_info(kmain,1)=J-2;
                orient(kmain,:)=new_orient_line;
                length_info(kmain,1)=norm(X(J,:)-X(J-2,:));
                age(kmain,1)=1; 
            end
            
            
            % Neighbor information for new curve
            if(change == 0)
                neigh(J-2,1) = -2;
                neigh(J-2,2) = J-1;
                neigh(J-1,1) = J-2;
                neigh(J-1,2) = J;
                neigh(J,1)   = J-1;
                neigh(J,2)   = -3;
                
            else
                neigh(J,1)   = -2;
                neigh(J,2)   = J-1;
                neigh(J-1,1) = J;
                neigh(J-1,2) = J-2;
                neigh(J-2,1) = J-1;
                neigh(J-2,2) = -3;
            end
            % Set new lines in Lambda
            n_triple = size(Lambda,1);
            Lambda(n_triple+1,:) = [i11,i12,J-2,k11,k12,kmain];
            Lambda(n_triple+2,:) = [i21,i22,J,  k21,k22,kmain];
        end
        %% Delete index_info, age and length_info entry of additional deleted subcurves
        % (2 subcurves at triple junction --> 1 subcurve, in add_delete is the
        % main- and sub-curve number of the sub-curve to be deleted, only the
        % entries in index_info and length_info, not the nodes!
        mark_list_new = mark_list;
        for i=1:add_delete_count
            k=add_delete(1,i);
            l0=add_delete(2,i);
            nr_curves(k,1)=nr_curves(k,1)-1;
            index_info(k,l0)=0; % set to 0, sort row later
            age(k,l0)=0;        % set to 0, sort row later
            length_info(k,l0)=-1; % set to -1, sort row later
            
            % Set entry in mark_list to -1
            row=find(mark_list((mark_count+1):mark_list_count,1)==k);
            Nrow = size(row,1);
            for ii = 1:Nrow
                row(ii)
                if(mark_list(row(ii)+mark_count,2)==l0)
                    mark_list_new(row(ii)+mark_count,:)=[-1,-1];
                end
            end
            for ii=1:mark_list_count
                if(mark_list(ii,1)==k && mark_list(ii,2)>l0)
                    mark_list_new(ii,2)=mark_list_new(ii,2)-1;
                end
            end
            
        end
        
        
        % sort row, shift 0 (-1) entries to the end
        for k=1:size(index_info,1)
            index_info(k,:)=sort_entries(index_info(k,:), 0);
            age(k,:)=sort_entries(age(k,:), 0); 
            length_info(k,:)=sort_entries(length_info(k,:), -1);
        end
        
        
        %% Delete now curve
        % help quantities
        J=size(X,1);
        X_new = X;
        neigh_new = neigh;
        index_info_new = index_info;
        nr_curves_new = nr_curves;
        age_new = age; 
        length_info_new = length_info;
        orient_new = orient;
        Lambda_new = Lambda;
        if(dim==3)
            closest_simp_new = closest_simp;
        end
        
        N0 = size(index_info,1);
        % Reduce numbers in neigh and index_info and Lambda where number > number of deleted point(s)
        for k=1:nr_points_deleted
            for i=1:J
                if(neigh(i,1)>info(k))
                    neigh_new(i,1)=neigh_new(i,1)-1;
                end
                if(neigh(i,2)>info(k))
                    neigh_new(i,2)=neigh_new(i,2)-1;
                end
            end
            for lll=1:N0
                for l0=1:nr_curves(lll,1)
                    if(index_info(lll,l0)>info(k))
                        index_info_new(lll,l0)=index_info_new(lll,l0)-1;
                    end
                end
            end
            for i=1:size(Lambda,1)
                for j=1:3
                    if(Lambda(i,j)>info(k))
                        Lambda_new(i,j)=Lambda_new(i,j)-1;
                    end
                end
            end
        end
        
        
        % Delete rows in X_new and neigh_new
        if(nr_points_deleted > 0)
            X_new(info(1:nr_points_deleted),:)=[];
            neigh_new(info(1:nr_points_deleted),:)=[];
            if(dim==3)
                closest_simp_new(info(1:nr_points_deleted),:) = [];
            end
                
        end
        % Delete rows in Lambda
        if(size(Lambda_new,1)>0 && row_delete_Lambda_count > 0)
            Lambda_new(row_delete_Lambda,:)=[];
        end
        
        % Delete index_info, length_info, nr_curves and age entry
        index_info_new(ll,l)  = 0;
        index_info_new(ll,:)  = sort_entries(index_info_new(ll,:),0);
        age_new(ll,l)             = 0; 
        age_new(ll,:)             = sort_entries(age_new(ll,:), 0); 
        length_info_new(ll,l) = -1;
        length_info_new(ll,:) = sort_entries(length_info_new(ll,:),-1);
        nr_curves_new(ll,1)   = nr_curves(ll,1) - 1;
        
        % Reduce number in mark_list
        for ii=1:mark_list_count
            if(mark_list(ii,1)==ll && mark_list(ii,2)>l)
                mark_list_new(ii,2)=mark_list_new(ii,2)-1;
            end
        end
        mark_list = mark_list_new;
        
        %% Update Gamma 
        X = X_new;
        neigh = neigh_new;
        index_info = index_info_new;
        age = age_new; 
        length_info = length_info_new;
        nr_curves = nr_curves_new;
        orient = orient_new;
        Lambda = Lambda_new;
        if(dim==3)
            closest_simp = closest_simp_new;
        end
    end
    
end

%% Delete entry if index_info(row,1)=0  i.e. main curve has to be deleted! 
mark0 = zeros(size(index_info_new,1),1);
count = 0;
for k=1:size(index_info,1);
    if(index_info_new(k,1)==0)
        count = count + 1;
        mark0(count)=k;
    end
end


index_info_new(mark0(1:count),:)=[];
length_info_new(mark0(1:count),:)=[];
nr_curves_new(mark0(1:count),:)=[];
age_new(mark0(1:count),:)=[]; 
orient_new(mark0(1:count),:)=[];

% Adapt Lambda(:,4:6) where main curve indices are stored!
Lambda_new2 = Lambda_new; 
for k=1:count
    kdel = mark0(k); 
    for i=1:size(Lambda_new,1)
        for j=4:6
            if(Lambda_new(i,j)>kdel)
                Lambda_new2(i,j)=Lambda_new2(i,j)-1; % reduce index by 1
            end
        end
    end
end
Lambda_new = Lambda_new2; 
        


%% Update Gamma 
Gamma.X = X_new;
Gamma.neigh = neigh_new;
Gamma.index_info = index_info_new;
Gamma.length = length_info_new;
Gamma.nr_curves = nr_curves_new;
Gamma.orient = orient_new;
Gamma.Lambda = Lambda_new;
Gamma.X_old  = Gamma.X;
Gamma.age    = age_new; 
Gamma.mark   = zeros(size(Gamma.index_info));

if(dim==3)
    Gamma.closest_simp = closest_simp; 
end

%% Write output
Gamma_new = Gamma;

end
