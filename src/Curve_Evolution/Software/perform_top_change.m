function [Gamma_new,Prev_new] = perform_top_change(split,merge,triple,boundary,Gamma,Image,Surface,config,Prev)
%% Initializations / Preparations 
J = size(Gamma.X,1); 
X = Gamma.X; 
neigh = Gamma.neigh; 
index_info = Gamma.index_info; 
nr_curves = Gamma.nr_curves; 
age = Gamma.age; 
orient = Gamma.orient; 
Lambda = Gamma.Lambda; 
I = Image.data; 
dim = size(Gamma.X,2); 
if(dim==3)
    closest_simp = Gamma.closest_simp;
end

% Help matrix needed if curve index is overwritten/changed due to previous
% top changes, matrix contains indexes of main and sub curve which the
% node belongs to
info_curve = zeros(J,2);   

% Matrix with information about previous top changes
Prev_new = Prev;

%% Boundary points - only 2d
nr_points = size(boundary,1); 
for count = 1:nr_points
        
    % Index of new boundary point
    i = boundary(count,1); 
        
    % Is this a point where a top change has happened in the last time
    % step? 
    p = X(i,:);
    Prev_p = in_Prev(p,Prev,config); 
    

    if(~Prev_p)
        fprintf('Topological changes occured! Boundary contact. \n');

        
        % Main and sub curve index
        if(info_curve(i,1)>0)
            ll = info_curve(i,1);
            l  = info_curve(i,2); 
        else
            ll = boundary(count,2);
            l  = boundary(count,3); 
        end
        
        % Boundary index
        boundary_index = boundary(count,4); 
        
        % Increase nr of nodes by 1 
        J=J+1; 
        
        % Change neighbor info
        i_m = neigh(i,1); 
        neigh(i,1) = -1;
        neigh(i_m,2) = J;
        neigh(J,1) = i_m;
        neigh(J,2) = -1; 
        
        % Project node to boundary
        switch boundary_index
            case 1 
                X(i,2)=0.5;
            case 2
                X(i,1)=size(I,2)+0.5;
            case 3
                X(i,2)=size(I,1)+0.5;
            case 4
                X(i,1)=0.5;
        end
        
        % Duplicate X(i,:)
        X(J,:)=X(i,:); 
        
        % Change index_info(ll,l)
        if(index_info(ll,l)~=i)
            if(neigh(index_info(ll,l),1)>0)  % Former closed curve
                index_info(ll,l)=i;
            else              % Former open curve, increase curve number
                nr_curves(ll,1) = nr_curves(ll,1) + 1;
                index_info(ll,nr_curves(ll,1))=i;
                age(ll,nr_curves(ll,1))=age(ll,l); 
                
                % Update info_curve 
                i0=i;
                count0 = 1;
                while(i0~= i || count0==1)
                    info_curve(i0,:)=[ll,nr_curves(ll,1)];
                    
                    i0 = neigh(i0,2); 
                    if(i0<0)
                        i0=i;
                    end
                    count0 = count0+1;
                end
            end
        end
        
        % Update Prev_new 
        n_Prev = size(Prev_new,1); 
        Prev_new(n_Prev+1,:)=p;
    end
        
end


%% Splitting points
nr_points = size(split,1); 

for count = 1:nr_points
                        
    % Index of nodes where splitting occurs
    i = split(count,1); 
    j = split(count,2); 
    
    
    % Is this a point where a top change has happened in the last time
    % step? 
    p = X(i,:); 
    Prev_p = in_Prev(p,Prev,config); 

    if(~Prev_p)
        
        % Main and sub curve index
        if(info_curve(i,1)>0)
            ll = info_curve(i,1);
            l  = info_curve(i,2); 
        else
            ll = split(count,3); 
            l  = split(count,4); 
        end
        
        % Check if curve indices are still the same, perform splitting only
        % in this case 
        if(info_curve(i,1)==info_curve(j,1) && info_curve(i,2)==info_curve(j,2))
            fprintf('Topological changes occured! Splitting of a curve. \n');
        
            % Insert two new points
            J=J+2; 
            X(J-1,:)=(X(i,:)+X(neigh(j,2),:))/2;
            X(J,:)  =(X(j,:)+X(neigh(i,2),:))/2; 
            
            
            % project to surface if dimension_flag = 3
            if(dim == 3)
                [closest_simp(J-1,1),X(J-1,:)] = get_closest_simp_use_parents(X(J-1,:),closest_simp(i,1),Surface);
                [closest_simp(J  ,1),X(J  ,:)] = get_closest_simp_use_parents(X(J  ,:),closest_simp(j,1),Surface);
                    
                               

            end

       
            % Change neighbor information (as in Paus and Benes)
            i_p = neigh(i,2); 
            j_p = neigh(j,2); 
            
            % new line j - J - i_p,   i - (J-1) - j_p
            neigh(j,2)  =J;
            neigh(J,1)  =j;
            neigh(J,2)  =i_p;
            neigh(i_p,1)=J;
            
            neigh(i,2)  =J-1;
            neigh(J-1,1)=i;
            neigh(J-1,2)=j_p;
            neigh(j_p,1)=J-1; 
                        
            % Adapt curve number
            nr_curves(ll,1)=nr_curves(ll,1)+1;

            index_info(ll,nr_curves(ll,1))=i; 

            age(ll,nr_curves(ll,1))=age(ll,l); 
            
            % Update info curve
            i0=index_info(ll,nr_curves(ll,1));
            count0 = 1;
            while(i0~= index_info(ll,nr_curves(ll,1)) || count0==1)
                info_curve(i0,:)=[ll,nr_curves(ll,1)];
                    
                i0 = neigh(i0,2); 
                if(i0<0)
                    i0=index_info(ll,nr_curves(ll,1));
                end
                count0 = count0+1;
            end

        
            % Update Prev_new
            n_Prev = size(Prev_new,1); 
            Prev_new(n_Prev+1,:)=p;

  
        end
    end
end

        
        
%% Merging points
nr_points = size(merge,1); 

for count = 1:nr_points
                        
    % Index of nodes where splitting occurs
    i = merge(count,1); 
    j = merge(count,4); 
     
        
    % Is this a point where a top change has happened in the last time
    % step? 
    p = X(i,:); 
    Prev_p = in_Prev(p,Prev,config); 

    if(~Prev_p)
        
        % Main and sub curve index
        if(info_curve(i,1)>0)
            ll = info_curve(i,1);
            l  = info_curve(i,2); 
        else
            ll = merge(count,2); 
            l  = merge(count,3); 
        end
        if(info_curve(j,1)>0)
            kk = info_curve(j,1);
            k  = info_curve(j,2); 
        else
            kk = merge(count,5); 
            k  = merge(count,6); 
        end
        
        
        fprintf('Topological changes occured! Merging of curves. \n');
        
        % Insert two new points
        J=J+2; 
        X(J-1,:)=(X(i,:)+X(neigh(j,2),:))/2;
        X(J,:)  =(X(j,:)+X(neigh(i,2),:))/2; 
        
        % project to surface if dimension_flag = 3
        if(dim == 3)
            [closest_simp(J-1,1),X(J-1,:)] = get_closest_simp_use_parents(X(J-1,:),closest_simp(i,1),Surface);
            [closest_simp(J  ,1),X(J  ,:)] = get_closest_simp_use_parents(X(J  ,:),closest_simp(j,1),Surface);
            
        end


        % Change neighbor information (as in Paus and Benes)
        i_p = neigh(i,2); 
        j_p = neigh(j,2); 

        % new line j - J - i_p,   i - (J-1) - j_p
        neigh(j,2)  =J;
        neigh(J,1)  =j;
        neigh(J,2)  =i_p;
        neigh(i_p,1)=J;

        neigh(i,2)  =J-1;
        neigh(J-1,1)=i;
        neigh(J-1,2)=j_p;
        neigh(j_p,1)=J-1; 

        
        % Delete entry in index_info for one of the curve
        % - If curve corresponding to index_info(kk,k) is open and curve
        %   corresponding to index_info(ll,l) is closed, delete
        %   index_info(kk,k)
        % - If curve corresponding to index_info(kk,k) is closed, delete
        %   index_info(kk,k), (curve ll,l arbitrary)
        % - Delete none of the curves if both are open. 
        
        first_l = index_info(ll,l); 
        first_k = index_info(kk,k); 
        
        if(neigh(first_l,1)>0 && neigh(first_k,1)<0)
            % Delete index_info entry (ll,l), i.e. nodes of (ll,l) belong
            % to (kk,k) in the following
            index_info(ll,l) = 0;
            % Delete age entry (ll,l)
            age(ll,l) = 0; 
            
            % Update info curve
            i0=index_info(kk,k);
            count0 = 1;
            while(i0~= index_info(kk,k) || count0==1)
                info_curve(i0,:)=[kk,k];
                i0 = neigh(i0,2); 
                if(i0<0)
                    i0=index_info(kk,k);
                end
                count0 = count0+1;
            end
        else
            if(neigh(first_k,1)>0)
                % Delete index_info entry (kk,k), i.e. nodes of (kk,k)
                % belong to (ll,l) in the following
                index_info(kk,k) = 0;
                % Delete age entry (kk,k)
                age(kk,k) = 0; 
                
                % Update info curve
                i0=index_info(ll,l);
                count0 = 1;

                while(i0~= index_info(ll,l) || count0==1)

                    info_curve(i0,:)=[ll,l];
                    i0 = neigh(i0,2); 
                    if(i0<0)
                        i0=index_info(ll,l);
                    end
                    count0 = count0+1;
                end
            else
                % Delete none of the curves index_info entries, but update
                % info curve
                i0=index_info(kk,k);
                count0 = 1;
                while(i0~= index_info(kk,k) || count0==1)
                    info_curve(i0,:)=[kk,k];
                    i0 = neigh(i0,2); 
                    if(i0<0)
                        i0=index_info(kk,k);
                    end
                    count0 = count0+1;
                end
                i0=index_info(ll,l);
                count0 = 1;
                while(i0~= index_info(ll,l) || count0==1)
                    info_curve(i0,:)=[ll,l];
                    i0 = neigh(i0,2); 
                    if(i0<0)
                        i0=index_info(ll,l);
                    end
                    count0 = count0+1;
                end
            end
               
        end
        % Note: index_info contains zero entries, the remaining entries are 
        % shifted at the end of the function (affects index_info and
        % orient)
        
        
        % Update Prev_new
        n_Prev = size(Prev_new,1); 
        Prev_new(n_Prev+1,:)=p;

    end
end

%% Triple junctions 
nr_points = size(triple,1); 

for count = 1:nr_points
                        
    % Index of nodes where splitting occurs
    i = triple(count,1); 
    j = triple(count,4); 
        
    % Is this a point where a top change has happened in the last time
    % step? 

    p = X(i,:); 
    Prev_p = in_Prev(p,Prev,config); 

    if(~Prev_p)
        
        % Main and sub curve index
        if(info_curve(i,1)>0)
            ll = info_curve(i,1);
            l  = info_curve(i,2); 
        else
            ll = triple(count,2); 
            l  = triple(count,3); 
        end
        if(info_curve(j,1)>0)
            kk = info_curve(j,1);
            k  = info_curve(j,2); 
        else
            kk = triple(count,5); 
            k  = triple(count,6); 
        end
        
        % Get phase / orientation info
        orient_l = zeros(1,2);
        orient_l(1) = orient(ll,1);
        orient_l(2) = orient(ll,2); 
        orient_k = zeros(1,2);
        orient_k(1) = orient(kk,1);
        orient_k(2) = orient(kk,2); 
        
        
        fprintf('Topological changes occured! Creation of triple junctions. \n');
        
        % Change neighbor information, perform triple junctions creation       
        case0 = 0; 
        if(orient_l(1)==orient_k(1))
            case0 = 1;
        else
            if(orient_l(1)==orient_k(2))
                case0 = 2; 
            else
                if(orient_l(2)==orient_k(1))
                    case0 = 3;
                else
                    if(orient_l(2)==orient_k(2))
                        case0 = 4;
                    else 
                        fprintf('Error: no matching regarding orient(i0,1/2) and orient(j0,1/2)\n');
                        pause
                    end
                end
            end
        end
        
        switch case0
            case 1
                % orient_l(1)==orient_k(1)
                h1 = norm(X(neigh(i,2),:)-X(neigh(j,1),:)); 
                h2 = norm(X(neigh(i,1),:)-X(neigh(j,2),:));
                if(h1<h2)
                    i1 = neigh(i,2); 
                    j1 = neigh(j,1); 
                    info_triple = [i1,i,j,j1,orient_k(2),orient_l(2)];
                else
                    i1 = neigh(i,1);
                    j1 = neigh(j,2); 
                    info_triple = [i,i1,j1,j,orient_l(2),orient_k(2)];
                end
        
            case 2
                % orient_l(1)==orient_k(2)
                h1 = norm(X(neigh(i,2),:)-X(neigh(j,2),:));
                h2 = norm(X(neigh(i,1),:)-X(neigh(j,1),:)); 
                if(h1<h2)
                    i1 = neigh(i,2); 
                    j1 = neigh(j,2); 
                    info_triple = [i1,i,j1,j,orient_k(1),orient_l(2)];
                else
                    i1 = neigh(i,1);
                    j1 = neigh(j,1); 
                    info_triple = [i,i1,j,j1,orient_l(2),orient_k(1)];
                end
            case 3
                % orient_l(2)==orient_k(1)
                h1 = norm(X(neigh(i,1),:)-X(neigh(j,1),:));
                h2 = norm(X(neigh(i,2),:)-X(neigh(j,2),:)); 
                if(h1<h2)
                    i1 = neigh(i,1); 
                    j1 = neigh(j,1); 
                    info_triple = [i,i1,j,j1,orient_k(2),orient_l(1)];
                else
                    i1 = neigh(i,2);
                    j1 = neigh(j,2); 
                    info_triple = [i1,i,j1,j,orient_l(1),orient_k(2)];
                end
            case 4
                % orient_l(2)==orient_k(2)
                h1 = norm(X(neigh(i,1),:)-X(neigh(j,2),:));
                h2 = norm(X(neigh(i,2),:)-X(neigh(j,1),:)); 
                if(h1<h2)
                    i1 = neigh(i,1); 
                    j1 = neigh(j,2); 
                    info_triple = [i,i1,j1,j,orient_k(1),orient_l(1)];
                else
                    i1 = neigh(i,2);
                    j1 = neigh(j,1); 
                    info_triple = [i1,i,j,j1,orient_l(1),orient_k(1)];
                end
        end
        % Change neighbor information of existing points
        neigh_old = neigh; 
        neigh(info_triple(1),1) = -2; 
        neigh(info_triple(2),2) = -3; 
        neigh(info_triple(3),1) = -2;
        neigh(info_triple(4),2) = -3;
        
        % Project nodes to center point
        p_Lambda1 = (X(i,:)+X(j,:))/2;
        p_Lambda2 = (X(i1,:)+X(j1,:))/2; 
        
        % project to surface if dimension_flag = 3
        if(dim == 3)
            [i_closest_simp1,p_Lambda1] = get_closest_simp_use_parents(p_Lambda1,closest_simp(i,1),Surface);
            [i_closest_simp2,p_Lambda2] = get_closest_simp_use_parents(p_Lambda2,closest_simp(i1,1),Surface);
            
            % Set new entrie closest_simp
            closest_simp(i ,1) = i_closest_simp1;
            closest_simp(i1,1) = i_closest_simp2;
            closest_simp(j ,1) = i_closest_simp1;
            closest_simp(j1,1) = i_closest_simp2;

        end
        
        
        % Set new points
        X(i,:)  = p_Lambda1; 
        X(i1,:) = p_Lambda2; 
        X(j,:)  = p_Lambda1;
        X(j1,:) = p_Lambda2; 
        
        
        
        
        % Create 3 new points 
        J=J+3; 
        X(J-2,:) = p_Lambda1;
        X(J-1,:) = (p_Lambda1 + p_Lambda2)/2;
        X(J,:)   = p_Lambda2;
        
        if(dim==3)
            closest_simp(J-2,1) = i_closest_simp1;
            closest_simp(J  ,1) = i_closest_simp2; 
            
            [closest_simp(J-1,1),X(J-1,:)] = get_closest_simp_use_parents(X(J-1,:),i_closest_simp1,Surface); 
        end

        % Neighbor information for the new points
        neigh(J-2,1) = -2; 
        neigh(J-2,2) = J-1; 
        neigh(J-1,1) = J-2;
        neigh(J-1,2) = J;
        neigh(J,1)   = J-1; 
        neigh(J,2)   = -3; 
            
            
        
        % Update 1st part of Lambda
        n_triple = size(Lambda,1); 
        n_triple = n_triple + 2; 
        
        Lambda(n_triple-1,1:3) = [i,j,J-2];      % Lambda_1
        Lambda(n_triple,1:3)   = [i1,j1,J];      % Lambda_2
        
           
       
        % Update index_info of existing curves
        l0 = index_info(ll,l); 
        k0 = index_info(kk,k); 
        
        
        if(neigh_old(l0,1)>0)  % closed curve
            index_info(ll,l) = info_triple(1);
        else                   % open curve, create new curve
            nr_curves(ll,1) = nr_curves(ll,1) + 1; 
            index_info(ll,nr_curves(ll,1)) = info_triple(1); 
            age(ll,nr_curves(ll,1)) = age(ll,l); 
            % Update info_curve
            i0=index_info(ll,nr_curves(ll,1));
            count0 = 1;
            while(i0~= index_info(ll,nr_curves(ll,1)) || count0==1)
                info_curve(i0,:)=[ll,nr_curves(ll,1)];

                i0 = neigh(i0,2); 
                if(i0<0)
                    i0=index_info(ll,nr_curves(ll,1));
                end
                count0 = count0+1;
            end
        end
        
        if(neigh_old(k0,1)>0)  % closed curve
            index_info(kk,k) = info_triple(3);
        else                   % open curve, create new curve
            nr_curves(kk,1) = nr_curves(kk,1) + 1; 
            index_info(kk,nr_curves(kk,1)) = info_triple(3); 
            age(kk,nr_curves(kk,1)) = age(kk,k);
           
            % Update info_curve
            i0=index_info(kk,nr_curves(kk,1));
            count0 = 1;
            while(i0~= index_info(kk,nr_curves(kk,1)) || count0==1)
                info_curve(i0,:)=[kk,nr_curves(kk,1)];

                i0 = neigh(i0,2); 
                if(i0<0)
                    i0=index_info(kk,nr_curves(kk,1));
                end
                count0 = count0+1;
            end
        end
        
        % Create new curve, new entry in index_info and orient
        new_orient_line = [info_triple(5), info_triple(6)];

        N0 = size(orient,1); 
        count0 = 1; 
        found = 0; 
        change = 0; 
        while(count0 <= N0 && found == 0)  %search for main curve index (in case that interface already exists)
            if(orient(count0,1)==new_orient_line(1) && orient(count0,2)==new_orient_line(2))
                % Interface already exists
                found = 1; 
                kmain = count0; 
                nr_curves(kmain,1) = nr_curves(kmain,1)+1;
                index_info(kmain,nr_curves(kmain,1)) = J-2; 
                age(kmain,nr_curves(kmain,1)) = 1; 
            else 
                if(orient(count0,1)==new_orient_line(2) && orient(count0,2)==new_orient_line(1))
                    found = 1; 
                    change = 1; % change orientation of the new curve from J-2,J-1,J to J,J-1,J
                    kmain = count0;
                    nr_curves(kmain,1)=nr_curves(kmain,1)+1;
                    index_info(kmain,nr_curves(kmain,1)) = J; % !!!!
                    age(kmain,nr_curves(kmain,1)) = 1; 
                end
            end
            count0 = count0 + 1; 
        end
        if(found==0)
            kmain = N0+1; 
            nr_curves(kmain,1) = 1; 
            index_info(kmain,nr_curves(kmain,1)) = J-2; 
            age(kmain,nr_curves(kmain,1))=1; 
            orient(kmain,1) = new_orient_line(1); 
            orient(kmain,2) = new_orient_line(2); 
        end
        
        % Change neighbor information (= orientation) if change==1
        % Neighbor information for the new points
        if(change==1)
            neigh(J,1)   = -2;
            neigh(J,2)   = J-1;
            neigh(J-1,1) = J;
            neigh(J-1,2) = J-2;
            neigh(J-2,1) = J-1;
            neigh(J-2,2) = -3;
        end
        
        
        
        % Update 2nd part of Lambda
        Lambda(n_triple-1,4:6) = [ll,kk,kmain];      % Lambda_1
        Lambda(n_triple,4:6)   = [ll,kk,kmain];      % Lambda_2
        
        % Update Prev_new
        n_Prev = size(Prev_new,1); 
        Prev_new(n_Prev+1,:)=p;
    end
end
            
        
%% Shift entries in index_info, age and orient, decrease nr_curves
N0 = size(index_info,1); 
for ll=1:N0
    N1 = nr_curves(ll,1);
    count_ll = 0; 
    collect = zeros(0,1); 
    for l = 1:N1
        if(index_info(ll,l)==0)
            count_ll            = count_ll + 1;
            collect(count_ll,1) = l; 
        end
    end
    for i = 1:count_ll
        l = collect(i)-i+1; 
        index_info(ll,l:(N1-1)) = index_info(ll,(l+1):N1); 
        index_info(ll,N1)       = 0;
        age(ll,l:(N1-1))        = age(ll,(l+1):N1); 
        age(ll,N1)              = 0; 
    end
        
    nr_curves(ll,1)=nr_curves(ll,1)-count_ll;    
end
     
%% Set output
Gamma_new            = Gamma; 
Gamma_new.X          = X; 
Gamma_new.neigh      = neigh; 
Gamma_new.index_info = index_info;
Gamma_new.nr_curves  = nr_curves;
Gamma_new.age        = age; 
Gamma_new.orient     = orient; 
Gamma_new.Lambda     = Lambda;
Gamma_new.X_old      = Gamma_new.X; 
Gamma_new.mark       = zeros(size(Gamma_new.index_info)); 

if(dim==3)
    Gamma_new.closest_simp = closest_simp;
end
Gamma_new = calc_normal_field(Gamma_new,Surface);

end


function out = in_Prev(p,Prev,config)
n = size(Prev,1); 

found = 0; 
j=1; 
while(j<=n && found==0)
    q = Prev(j,:); 
    if(norm(q-p)<config.top_check.a_Prev)
        found = 1; 
    end
    j=j+1;
end

out = found;
end
