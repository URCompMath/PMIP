function Gamma = get_initial_curve(config,initial_information,Surface)

%% Initial index_info
Nmain = size(initial_information.radius,3);
Nsub  = size(initial_information.radius,2);

Gamma.index_info = zeros(Nmain,Nsub);
Gamma.nr_curves = Nsub * ones(Nmain,1);
Gamma.age = zeros(Nmain,Nsub);

J0 = round(initial_information.J / (Nmain*Nsub));

for i=1:Nmain
    for j=1:Nsub
        Gamma.index_info(i,j) = (i-1)*Nsub*J0 + (j-1)*J0 + 1;
    end
end

%% Initial curve nodes
% Assume closed curves with center c and radius r at the beginning

X = zeros(initial_information.J,config.dimension);
c = initial_information.center;
r = initial_information.radius;

if(config.dimension==3)
    closest_simp = zeros(initial_information.J,1);
    i_parent = -1;
end

for j=1:Nmain
    for k=1:Nsub
        Gamma.age(j,k)=1;
        i0 = Gamma.index_info(j,k);
        if(k<Nsub)
            iN = Gamma.index_info(j,k+1)-1;
        else
            if(j<Nmain)
                iN=Gamma.index_info(j+1,1)-1;
            else
                iN=initial_information.J;
            end
        end
        J_curve = iN-i0+1;

        for i=i0:iN
            switch config.dimension
                case 2 % Curves in R^2
                    X(i,1) = c(1,k,j) + r(1,k,j)*cos(2*pi*(i-i0+1)/J_curve);
                    X(i,2) = c(2,k,j) + r(2,k,j)*sin(2*pi*(i-i0+1)/J_curve);
                    
                    
                case 3 % Curves on surfaces (nodes in 3D, but restricted on surface)
                    % Create initial curves before projection to surface
                    switch Surface.flag
                            
                        case 1           % Bunny
                            X(i,1) = c(1,k,j) + r(1,k,j)*cos(2*pi*(i-i0+1)/J_curve);
                            X(i,2) = c(2,k,j) + r(2,k,j)*sin(2*pi*(i-i0+1)/J_curve);
                            X(i,3) = c(3,k,j);
                            
                            v = [0,0,-1];
                            
                        case 5          % Earth, longwave radiation example
                            if(k==1 || k==3)
                                X(i,1) = c(1,k,j);
                                X(i,2) = c(2,k,j) + (-1)^((k-1)/2) * r(2,k,j)*cos(2*pi*(i-i0+1)/J_curve);
                                X(i,3) = c(3,k,j) + r(3,k,j)*sin(2*pi*(i-i0+1)/J_curve);
                            else
                                if(k==2 || k==4)
                                    X(i,1) = c(1,k,j) + (-1)^(k/2) * r(1,k,j)*cos(2*pi*(i-i0+1)/J_curve);
                                    X(i,2) = c(2,k,j);
                                    X(i,3) = c(3,k,j) + r(3,k,j)*sin(2*pi*(i-i0+1)/J_curve);
                                    
                                else
                                    if(k==6)
                                        factor = -1;
                                    else
                                        factor = 1;
                                    end
                                    X(i,1) = c(1,k,j) + r(1,k,j)*cos(2*pi*(i-i0+1)/J_curve);
                                    X(i,2) = c(2,k,j) + factor * r(2,k,j)*sin(2*pi*(i-i0+1)/J_curve);
                                    X(i,3) = c(3,k,j);
                                    
                                end
                            end
                            
                            
                            v = -c(:,k,j)';
                            
                        case 6   % Earth, net radiation example                            
                            if(abs(c(3,k,j)) ~= 1)
                                alpha = atan2(c(2,k,j),c(1,k,j));
                                e1 = [-sin(alpha); cos(alpha); 0];
                                e2 = cross(e1,-c(:,k,j));
                            else
                                if(c(3,k,j)==1)
                                    e1 = [0;1;0];
                                    e2 = [-1;0;0];
                                else
                                    e1 = [0;1;0];
                                    e2 = [1;0;0];
                                end
                            end
                            rr = r(1,k,j);
                            
                            X(i,:) = c(:,k,j) + rr * cos(2*pi*(i-i0+1)/J_curve) * e1 + rr * sin(2*pi*(i-i0+1)/J_curve) * e2;
                            
                            v = -c(:,k,j)';
                            
                        case 7 % Torus
                            cc = c(:,k,j); 
                            cc = cc/norm(cc); 
                            e2 = [0;0;1]; 
                            e1 = cross(e2,cc); 
                            
                            X(i,:) = c(:,k,j) + r(1,k,j) * cos(2*pi*(i-i0+1)/J_curve) * (-e1) + r(2,k,j) * sin(2*pi*(i-i0+1)/J_curve) * e2; 
                            v = -c(:,k,j)';

                        otherwise        % Faces
                            X(i,1) = c(1,k,j) + r(1,k,j)*cos(2*pi*arc_fcn((i-i0+1)/J_curve));
                            X(i,2) = c(2,k,j) + r(2,k,j)*sin(2*pi*arc_fcn((i-i0+1)/J_curve));
                            X(i,3) = c(3,k,j);
                            
                            v = [0,0,-1];
                            
                            
                    end
                    
                    % Project nodes to Surface
                    if(i==i0)
                        [X(i,:),closest_simp(i,1)] = project2surface_init(X(i,:),v,Surface);
                    else
                        [X(i,:),closest_simp(i,1)] = project2surface_parent(X(i,:),v,i_parent,Surface);
                    end
                    i_parent = closest_simp(i,1);
                    
            end
        end
    end
end

Gamma.X  = X;
Gamma.X_old = X;

if(config.dimension == 3)
    Gamma.closest_simp = closest_simp;
end

%% Triple Junctions
% Assume no triple junctions at the beginning (only closed curves)
Gamma.Lambda = zeros(0,6);


%% Orientation of the curves
Gamma.orient = zeros(Nmain,2);
for i=1:Nmain
    Gamma.orient(i,1)=Nmain+1;
    Gamma.orient(i,2)=i;
end


end


function S=arc_fcn(s)

a = 2/3;

if(s<=1/8)
    S=2*a*s;
else
    if(s<=3/8)
        S = a*(s-1/8) + 1/6;
    else
        if(s<=5/8)
            S = 2*a*(s-3/8) + 1/3;
        else
            if(s<=7/8)
                S = a*(s-5/8) + 2/3;
            else
                S = 2*a*(s-7/8) + 5/6;
            end
        end
    end
end
end




