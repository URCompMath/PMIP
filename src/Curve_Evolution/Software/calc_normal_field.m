function Gamma_new = calc_normal_field(Gamma,Surface)

dim = size(Gamma.X,2);

% Intialize nu
Gamma.nu = zeros(size(Gamma.X,1),dim);

if(dim==3)
    Gamma.x_s  = zeros(size(Gamma.X,1),dim);
    Gamma.nu_F = zeros(size(Gamma.X,1),dim);
end

% Calculate normal to the curve, in case of curves on surfaces also compute
% tangent vector and normal vector to the surface
switch dim
    case 2
        for i=1:size(Gamma.X,1)
            if(Gamma.neigh(i,1)>0 && Gamma.neigh(i,2)>0)
                h1 = norm(Gamma.X_old(i,:)-Gamma.X_old(Gamma.neigh(i,1),:));
                h2 = norm(Gamma.X_old(Gamma.neigh(i,2),:)-Gamma.X_old(i,:));
                Gamma.nu(i,:)=([0,-1;1,0]*(Gamma.X_old(Gamma.neigh(i,2),:)-Gamma.X_old(Gamma.neigh(i,1),:))')';
                Gamma.nu(i,:)=Gamma.nu(i,:)/(h1+h2);
            else
                if(Gamma.neigh(i,1)<0)
                    if(Gamma.neigh(i,2)<0)
                        i
                        Gamma.neigh(i,:)
                    end
                    Gamma.nu(i,:)=([0,-1;1,0]*(Gamma.X_old(Gamma.neigh(i,2),:)-Gamma.X_old(i,:))')';
                    Gamma.nu(i,:)=Gamma.nu(i,:)/norm(Gamma.nu(i,:));
                else
                    Gamma.nu(i,:)=([0,-1;1,0]*(Gamma.X_old(i,:)-Gamma.X_old(Gamma.neigh(i,1),:))')';
                    Gamma.nu(i,:)=Gamma.nu(i,:)/norm(Gamma.nu(i,:));
                end
            end
        end
    case 3
        for i=1:size(Gamma.X,1)
            % Calc tangent vector
            if(Gamma.neigh(i,1)>0 && Gamma.neigh(i,2)>0)
                Gamma.x_s(i,:)=Gamma.X_old(Gamma.neigh(i,2),:)-Gamma.X_old(Gamma.neigh(i,1),:);
            else
                if(Gamma.neigh(i,1)<0)
                    Gamma.x_s(i,:)=Gamma.X_old(Gamma.neigh(i,2),:)-Gamma.X_old(i,:);
                else
                    Gamma.x_s(i,:)=Gamma.X_old(i,:)-Gamma.X_old(Gamma.neigh(i,1),:);
                end
            end
            Gamma.x_s(i,:)=Gamma.x_s(i,:)/norm(Gamma.x_s(i,:));
            % Calc normal vector to the surface
            Gamma.nu_F(i,:) = calc_normal_to_surface(i,Gamma,Surface);
            Gamma.nu(i,:) = -cross(Gamma.x_s(i,:), Gamma.nu_F(i,:)); 
            Gamma.nu(i,:) = Gamma.nu(i,:)/norm(Gamma.nu(i,:)); 
        end
        
        
end

Gamma_new = Gamma; 

end

