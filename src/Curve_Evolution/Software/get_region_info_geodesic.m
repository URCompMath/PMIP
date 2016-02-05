function Omega_new = get_region_info_geodesic(Gamma, Surface, Omega, Image, config, kmax, initial)
% This fcn sets Omega.A_info, Omega.I_info, Omega.n_info

Nsimp = size(Surface.tri,1);
Nmain = size(Gamma.index_info,1);

% Color flag
if(size(Image.data,2)==3)
    color_flag = 1;
else  
    color_flag = 0;
end

% help array
info = zeros(Nsimp,1);

% 1st step: Set triangles where the nodes lie to Gamma.orient(k,2)
nsimp_local = 0;
simp_local = zeros(nsimp_local,10);


for k=1:Nmain
    for l=1:Gamma.nr_curves(k,1)
        
        
        j=1;
        i=Gamma.index_info(k,l);
        while( i ~= Gamma.index_info(k,l) || j==1)
            i_simp = Gamma.closest_simp(i,1); 
            
            if(~info(i_simp,1))
                
                p1 = Surface.nodes(Surface.tri(i_simp,1),:);
                p2 = Surface.nodes(Surface.tri(i_simp,2),:);
                p3 = Surface.nodes(Surface.tri(i_simp,3),:);
                
                p = (p1+p2+p3)/3;
                x = Gamma.X(i,:);
                v = Gamma.nu(i,:);
                
                % Get I0
                I0 = Image.data(i_simp,:); 
                if(color_flag)
                    switch config.method.color
                        case 2 % CB
                            I0 = my_rgb2cb(I0)'; 
                        case 3 % HSV
                            I0 = my_rgb2hsv(I0)'; 
                            % I0(1) in [0,2pi] --> convert to S^2
                            I0 = [cos(I0(1)), sin(I0(1)), I0(2), I0(3)];
                    end
                end
                            
                            
                
                
                % Set A_info, I_info, n_info
                if((p-x)*v'>0)
                    if(initial)
                        Omega.A_info(i_simp,1) = Gamma.orient(k,2);
                        Omega.I_info(Gamma.orient(k,2),:) = Omega.I_info(Gamma.orient(k,2),:) + I0;
                        Omega.n_info(Gamma.orient(k,2),1) = Omega.n_info(Gamma.orient(k,2),1) + 1;
                        region_new = Gamma.orient(k,2);
                    else
                        region_old = Omega.A_info(i_simp,1);
                        region_new = Gamma.orient(k,2);
                        if(region_old ~= region_new)
                            Omega.A_info(i_simp,1) = region_new;
                            Omega.I_info(region_new,:) = Omega.I_info(region_new,:) + I0;
                            Omega.I_info(region_old,:) = Omega.I_info(region_old,:) - I0;
                            Omega.n_info(region_new,1) = Omega.n_info(region_new,:) + 1;
                            Omega.n_info(region_old,1) = Omega.n_info(region_old,:) - 1;
                        end
                    end
                else
                    if(initial)
                        Omega.A_info(i_simp,1) = Gamma.orient(k,1);
                        Omega.I_info(Gamma.orient(k,1),:) = Omega.I_info(Gamma.orient(k,1),:) + I0;
                        Omega.n_info(Gamma.orient(k,1),1) = Omega.n_info(Gamma.orient(k,1),1) + 1;
                        region_new = Gamma.orient(k,1);
                    else
                        region_old = Omega.A_info(i_simp,1);
                        region_new = Gamma.orient(k,1);
                        if(region_old ~= region_new)
                            Omega.A_info(i_simp,1) = region_new;
                            Omega.I_info(region_new,:) = Omega.I_info(region_new,:) + I0;
                            Omega.I_info(region_old,:) = Omega.I_info(region_old,:) - I0;
                            Omega.n_info(region_new,1) = Omega.n_info(region_new,:) + 1;
                            Omega.n_info(region_old,1) = Omega.n_info(region_old,:) - 1;
                        end
                    end
                end
                
                nsimp_local = nsimp_local + 1;
                simp_local(nsimp_local,:) = [i_simp,x,v,k,region_new,l];
                
                % Set help array
                info(i_simp,1) = 1;
            end
                    
                    
                
            i = Gamma.neigh(i,2);
            if(i<0)
                i = Gamma.index_info(k,l);
            end
            j=j+1;
        end
    end
end


kk=1;

kheritage = 5; 

% Loops: 
while(nsimp_local > 0 && kk<=kmax)
    % Next step: Look at the neighbors of simp_local
    nsimp_local_new = 0;
    simp_local_new = zeros(nsimp_local_new,10);
    
    for i = 1:nsimp_local
        for j=1:3
            i_neigh = Surface.neigh(simp_local(i,1),j);
            
            if(i_neigh > -1)
                % Consider neighbor triangle only if it has not been
                % considered yet:
                if(info(i_neigh,1)==0)
                    % Get I0
                    I0 = Image.data(i_neigh,:);
                    if(color_flag)
                        switch config.method.color
                            case 2 % CB
                                I0 = my_rgb2cb(I0)';
                            case 3 % HSV
                                I0 = my_rgb2hsv(I0)';
                                % I0(1) in [0,2pi] --> convert to S^2
                                I0 = [cos(I0(1)), sin(I0(1)), I0(2), I0(3)];
                            case 4 % Lab
                                % to be implemented
                        end
                    end
                       
                    if(kk<kheritage)
                        p1 = Surface.nodes(Surface.tri(i_neigh,1),:);
                        p2 = Surface.nodes(Surface.tri(i_neigh,2),:);
                        p3 = Surface.nodes(Surface.tri(i_neigh,3),:);
                        
                        p = (p1+p2+p3)/3;
                        x = simp_local(i,2:4);
                        v = simp_local(i,5:7);
                        k = simp_local(i,8);
                        ksub = simp_local(i,10); 
                        region_heritage = simp_local(i,9);


                
                        % Set A_info, I_info, n_info
                        if((p-x)*v'>0)  % inside
                            
                            % Set A_info, I_info, n_info
                            if(initial)
                                region_new = Gamma.orient(k,2);
                                Omega.A_info(i_neigh,1) = region_new;
                                Omega.I_info(region_new,:) = Omega.I_info(region_new,:) + I0;
                                Omega.n_info(region_new,1) = Omega.n_info(region_new,1) + 1;
                            else
                                region_old = Omega.A_info(i_neigh,1);
                                region_new = Gamma.orient(k,2);
                                if(region_old ~= region_new)
                                    Omega.A_info(i_neigh,1) = region_new;
                                    Omega.I_info(region_new,:) = Omega.I_info(region_new,:) + I0;
                                    Omega.I_info(region_old,:) = Omega.I_info(region_old,:) - I0;
                                    Omega.n_info(region_new,1) = Omega.n_info(region_new,:) + 1;
                                    Omega.n_info(region_old,1) = Omega.n_info(region_old,:) - 1;
                                end
                                
                                
                            end
                            % Set help array
                            info(i_neigh,1)=1;
                            
                            % Update simp_local list
                            nsimp_local_new = nsimp_local_new + 1;
                            
                            [x,v,k] = get_dist_and_normal_p(p,Gamma,k,ksub);
                            simp_local_new(nsimp_local_new,:) = [i_neigh,x,v,k,region_new,ksub];
                       
                        else
                            % outside
                            % Set A_info, I_info, n_info
                            if(initial)
                                region_new = Gamma.orient(k,1);
                                Omega.A_info(i_neigh,1) = region_new;
                                Omega.I_info(region_new,:) = Omega.I_info(region_new,:) + I0;
                                Omega.n_info(region_new,1) = Omega.n_info(region_new,1) + 1;
                            else
                                region_old = Omega.A_info(i_neigh,1);
                                region_new = Gamma.orient(k,1);
                                if(region_old ~= region_new)
                                    Omega.A_info(i_neigh,1) = region_new;
                                    Omega.I_info(region_new,:) = Omega.I_info(region_new,:) + I0;
                                    Omega.I_info(region_old,:) = Omega.I_info(region_old,:) - I0;
                                    Omega.n_info(region_new,1) = Omega.n_info(region_new,:) + 1;
                                    Omega.n_info(region_old,1) = Omega.n_info(region_old,:) - 1;
                                end
                            end
                            % Set help array
                            info(i_neigh,1)=1;
                            
                            % Update simp_local list
                            nsimp_local_new = nsimp_local_new + 1;
                            
                            [x,v,k] = get_dist_and_normal_p(p,Gamma,k,ksub);
                            simp_local_new(nsimp_local_new,:) = [i_neigh,x,v,k,region_new,ksub];
                           
                            
                        end
                    else
                        x = simp_local(i,2:4);
                        v = simp_local(i,5:7);
                        k = simp_local(i,8);
                        region_heritage = simp_local(i,9);
                        
                        if(initial)
                            region_new = region_heritage; 
                            Omega.A_info(i_neigh,1) = region_new; 
                            Omega.I_info(region_new,:) = Omega.I_info(region_new,:) + I0;
                            Omega.n_info(region_new,1) = Omega.n_info(region_new,1) + 1;
                        else
                            region_old = Omega.A_info(i_neigh,1);
                            region_new = region_heritage; 
                            if(region_old ~= region_new)
                                Omega.A_info(i_neigh,1) = region_new;
                                Omega.I_info(region_new,:) = Omega.I_info(region_new,:) + I0;
                                Omega.I_info(region_old,:) = Omega.I_info(region_old,:) - I0;
                                Omega.n_info(region_new,1) = Omega.n_info(region_new,:) + 1;
                                Omega.n_info(region_old,1) = Omega.n_info(region_old,:) - 1;
                            end
                        end
                        
                        % Set help array
                        info(i_neigh,1)=1;
                        
                        % Update simp_local list
                        nsimp_local_new = nsimp_local_new + 1;
                        simp_local_new(nsimp_local_new,:) = [i_neigh,x,v,k,region_new,ksub];

                    end
                end
            end
        end
    end
    
    % Update for next step kk+1
    simp_local = simp_local_new;
    nsimp_local = nsimp_local_new;
    kk = kk + 1;
    
       
end

Omega_new = Omega;
end


function [q,n,k] = get_dist_and_normal_p(p,Gamma,kmain,ksub)


dist = 1e8;
n = [0 0 0];
q = [0 0 0];
k=1;

i = Gamma.index_info(kmain,ksub);
j = 1;

while(i ~= Gamma.index_info(kmain,ksub) || j==1)
    xi = Gamma.X(i,:);
    if(Gamma.neigh(i,2)>-1)
        xii = Gamma.X(Gamma.neigh(i,2),:);
        
        v = xii-xi;
        l = norm(v);
        v = v/l;
        factor = (p-xi)*(v') ;
        if(factor >=0 && factor <= l)
            Pp = xi + ((p-xi)*(v'))*v;
            info = 0;
        else
            if(factor < 0)
                Pp = xi;
                info = 1;
            else
                Pp = xii;
                info = 2;
            end
        end
        disti = norm(p - Pp);
        if(disti < dist)
            dist = disti;
            k = kmain;
            switch info
                case 0
                    n1 = Gamma.nu(i,:);
                    n1 = n1/norm(n1);
                    n2 = Gamma.nu(Gamma.neigh(i,2),:);
                    n2 = n2/norm(n2);
                    n = (n1+n2)/2 ;
                case 1
                    n = Gamma.nu(i,:);
                case 2
                    n = Gamma.nu(Gamma.neigh(i,2),:);
            end
            n = n/norm(n);
            q = Pp;
            
        end
    end
    i = Gamma.neigh(i,2);
    if(i<0)
        i=Gamma.index_info(kmain,ksub);
    end
    j= j+1;
end

end
