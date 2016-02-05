function Image_new = image_diffusion_and_denoising_2d_sparse(Image,Omega,config)

%% Preparations
Nx = Image.sizes(1);
Ny = Image.sizes(2);

lambda = 1/config.param_diffusion;

h=1; 

A_info = Omega.A_info;
nphases = size(Omega.coeffs,1);

Help = zeros(Ny+1,Nx+1,nphases);

for i=1:Nx
    for j=1:Ny
        k = A_info(j,i);
        Help(j  ,i  ,k)=1;
        Help(j+1,i  ,k)=1;
        Help(j+1,i+1,k)=1;
        Help(j  ,i+1,k)=1;
    end
end


I0 = double(Image.data);

% Get color_flag
if(size(size(Image.data),2)==3)
    color_flag = 1;
else
    color_flag = 0;
end


if(color_flag == 0)
    u = zeros((Ny+1)*(Nx+1),nphases); 
    for k=1:nphases
        % Get sparse system matrix and right hand side
        fprintf('Phase %d\n', k); 
        [S,b] = get_matrix_rhs(Nx,Ny,I0,A_info,k,Help(:,:,k),lambda,h,0);
        % Solve
        u(:,k) = S\b;
        
    end
    fprintf('Get I_approx\n'); 
    I_approx = zeros(Ny,Nx); 
    for ix=1:Nx
        for iy=1:Ny
            k = A_info(iy,ix); 
            % (jx,ix)
            i = (ix-1)*(Ny+1)+iy;
            u_1 = u(i,k); 
            
            % (jx+1,ix)
            i = (ix-1)*(Ny+1)+iy+1;
            u_2 = u(i,k); 
            
            % (jx,ix+1)
            i = ix*(Ny+1)+iy;
            u_3 = u(i,k); 
            
            % (jx+1,ix+1)
            i = ix*(Ny+1)+iy+1;
            u_4 = u(i,k); 
            
            I_approx(iy,ix) = (u_1+u_2+u_3+u_4)/4;
        end
    end
else
    u1 = zeros((Ny+1)*(Nx+1),nphases); 
    u2 = zeros((Ny+1)*(Nx+1),nphases); 
    u3 = zeros((Ny+1)*(Nx+1),nphases); 
    for k=1:nphases
        % Get sparse system matrix and right hand side
        fprintf('Phase %d\n', k); 
        [S,b] = get_matrix_rhs(Nx,Ny,I0,A_info,k,Help(:,:,k),lambda,h,1);
        % Solve
        u1(:,k) = S\(b(:,1));
        u2(:,k) = S\(b(:,2));
        u3(:,k) = S\(b(:,3));
        
    end
    fprintf('Get I_approx\n'); 
    I_approx = zeros(Ny,Nx,3); 
     for ix=1:Nx
         for iy=1:Ny
             k = A_info(iy,ix);
             % (jx,ix)
             i = (ix-1)*(Ny+1)+iy;
             u_1 = [u1(i,k);u2(i,k);u3(i,k)];
             
             % (jx+1,ix)
             i = (ix-1)*(Ny+1)+iy+1;
             u_2 = [u1(i,k);u2(i,k);u3(i,k)];
             
             % (jx,ix+1)
             i = ix*(Ny+1)+iy;
             u_3 = [u1(i,k);u2(i,k);u3(i,k)];
             
             % (jx+1,ix+1)
             i = ix*(Ny+1)+iy+1;
             u_4 = [u1(i,k);u2(i,k);u3(i,k)];
             
             I_approx(iy,ix,:) = reshape((u_1+u_2+u_3+u_4)/4,1,1,3); 
             
         end
     end
end


%% Generate output
Image_new = Image;
Image_new.data = uint8(I_approx); 

% Draw final and original image
figure
if(color_flag == 0)  % scalar image, need colormap
    image(I_approx);
    colormap(Image.cmap);
    axis image
    
else
    image(uint8(I_approx));       % vector valued image, RGB
    axis image
end
pause(0.1);


end

function [S,b] = get_matrix_rhs(Nx,Ny,I0,A_info,phase_index,Help,lambda,h,color_flag)

N0 = 5*(Nx+1)*(Ny+1);
row = zeros(N0,1);
col = zeros(N0,1);
data = zeros(N0,1);

if(color_flag == 0)
    b = zeros((Nx+1)*(Ny+1),1);
else
    b = zeros((Nx+1)*(Ny+1),3); 
end


count = 0;

for iy=0:Ny
    for ix=0:Nx
        if(Help(iy+1,ix+1)==1)
            % use N(nord),E(east),S(south),W(west),
            
            N  = [ix  ,iy+1];
            S  = [ix  ,iy-1];
            E  = [ix+1,iy  ];
            W  = [ix-1,iy  ];
            
            NE = [ix+1,iy+1];
            SE = [ix+1,iy-1];
            SW = [ix-1,iy-1];
            NW = [ix-1,iy+1];
            
            
            tophase = zeros(4,1);
            neigh = [N;E;S;W];
            
            tophase_corners = zeros(4,1);
            neigh_corners = [NE;SE;SW;NW];
            
            for i=1:4
                jx = neigh(i,1);
                jy = neigh(i,2);
                if(jx>=0 && jx<=Nx && jy>=0 && jy<=Ny)
                    if(Help(jy+1,jx+1)==1)
                        tophase(i,1)=1;
                    end
                end
                
                jx = neigh_corners(i,1);
                jy = neigh_corners(i,2);
                if(jx>=0 && jx<=Nx && jy>=0 && jy<=Ny)
                    if(Help(jy+1,jx+1)==1)
                        tophase_corners(i,1)=1;
                    end
                end
            end
            
            %% Get type
            
            if(tophase(1) && tophase(2) && tophase(3) && tophase(4))
                % Look at the corners to differ between type I and type IV
                
                if(tophase_corners(1)==0) % type IV
                    % Point NE does not belong to phase with phase_index
                    % N and E --> -1
                    % S and W --> -2
                    
                    % star = C,N,E,S,W
                    star = [6 + 1.5*lambda*h^2, -1, -1, -2, -2];
                    rhs = 1.5*lambda*h^2*get_image_at_grid(ix,iy,I0,A_info,phase_index,color_flag);
                    
                else
                    if(tophase_corners(2)==0) % type IV
                        % Point SE does not belong to phase with phase_index
                        % S and E --> -1
                        % N and W --> -2
                        
                        % star = C,N,E,S,W
                        star = [6 + 1.5*lambda*h^2, -2, -1, -1, -2];
                        rhs = 1.5*lambda*h^2*get_image_at_grid(ix,iy,I0,A_info,phase_index,color_flag);
                        
                    else
                        if(tophase_corners(3)==0) % type IV
                            % Point SW does not belong to phase with phase_index
                            % S and W --> -1
                            % N and E --> -2
                            
                            % star = C,N,E,S,W
                            star = [6 + 1.5*lambda*h^2, -2, -2, -1, -1];
                            rhs = 1.5*lambda*h^2*get_image_at_grid(ix,iy,I0,A_info,phase_index,color_flag);
                            
                            
                        else
                            if(tophase_corners(4)==0) % type IV
                                % Point NW does not belong to phase with phase_index
                                % N and W --> -1
                                % S and E --> -2
                                
                                % star = C,N,E,S,W
                                star = [6 + 1.5*lambda*h^2, -1, -2, -2, -1];
                                rhs = 1.5*lambda*h^2*get_image_at_grid(ix,iy,I0,A_info,phase_index,color_flag);
                                
                            else % type I
                                % star = C,N,E,S,W
                                star = [8+2*lambda*h^2, -2, -2, -2, -2];
                                rhs  = 2*lambda*h^2 * get_image_at_grid(ix,iy,I0,A_info,phase_index,color_flag);
                            end
                        end
                    end
                end
            else
                if(tophase(1) && tophase(2) && tophase(3))  % type II, W missing
                    % Point W does not belong to phase
                    % N and S --> -1
                    % E       --> -2
                    % W       -->  0
                    
                    % star = C,N,E,S,W
                    star = [4+lambda*h^2, -1, -2, -1, 0];
                    rhs  = lambda*h^2 * get_image_at_grid(ix,iy,I0,A_info,phase_index,color_flag);
                    
                    
                else
                    if(tophase(1) && tophase(2) && tophase(4))  % type II, S missing
                        % Point S does not belong to phase
                        % E and W --> -1
                        % N       --> -2
                        % S       -->  0
                        
                        % star = C,N,E,S,W
                        star = [4+lambda*h^2, -2, -1, 0, -1];
                        rhs  = lambda*h^2 * get_image_at_grid(ix,iy,I0,A_info,phase_index,color_flag);
                        
                        
                    else
                        if(tophase(1) && tophase(3) && tophase(4)) % type II, E missing
                            % Point E does not belong to phase
                            % N and S --> -1
                            % W       --> -2
                            % E       -->  0
                            
                            % star = C,N,E,S,W
                            star = [4+lambda*h^2, -1, 0, -1, -2];
                            rhs  = lambda*h^2 * get_image_at_grid(ix,iy,I0,A_info,phase_index,color_flag);
                            
                            
                            
                        else
                            if(tophase(2) && tophase(3) && tophase(4)) % type II, N missing
                                % Point N does not belong to phase
                                % E and W --> -1
                                % S       --> -2
                                % N       -->  0
                                
                                % star = C,N,E,S,W
                                star = [4+lambda*h^2, 0, -1, -2, -1];
                                rhs  = lambda*h^2 * get_image_at_grid(ix,iy,I0,A_info,phase_index,color_flag);
                                
                                
                            else
                                if(tophase(1) && tophase(2))  %% type III, S and W missing
                                    % Points S and W do not belong to phase
                                    % N and E --> -1
                                    % S and W -->  0
                                    
                                    % star = C,N,E,S,W
                                    star = [2+0.5*lambda*h^2, -1, -1, 0, 0];
                                    rhs  = 0.5*lambda*h^2 * get_image_at_grid(ix,iy,I0,A_info,phase_index,color_flag);
                                    
                                else
                                    if(tophase(1) && tophase(4))  %% type III, S and E missing
                                        % Points S and E do not belong to phase
                                        % N and W --> -1
                                        % S and E -->  0
                                        
                                        % star = C,N,E,S,W
                                        star = [2+0.5*lambda*h^2, -1, 0, 0, -1];
                                        rhs  = 0.5*lambda*h^2 * get_image_at_grid(ix,iy,I0,A_info,phase_index,color_flag);
                                        
                                        
                                    else
                                        if(tophase(2) && tophase(3))  %% type III, N and W are missing
                                            % Points N and W do not belong to phase
                                            % S and E --> -1
                                            % N and W -->  0
                                            
                                            % star = C,N,E,S,W
                                            star = [2+0.5*lambda*h^2, 0, -1, -1, 0];
                                            rhs  = 0.5*lambda*h^2 * get_image_at_grid(ix,iy,I0,A_info,phase_index,color_flag);
                                            
                                            
                                            
                                        else
                                            if(tophase(3) && tophase(4))  %% type III, N and E missing
                                                % Points N and E do not belong to phase
                                                % S and W --> -1
                                                % N and E -->  0
                                                
                                                % star = C,N,E,S,W
                                                star = [2+0.5*lambda*h^2, 0, 0, -1, -1];
                                                rhs  = 0.5*lambda*h^2 * get_image_at_grid(ix,iy,I0,A_info,phase_index,color_flag);
                                                
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                
                
            end
            
        
        
            r  = ix*(Ny+1)+iy+1; % row (center)
            cN = ix*(Ny+1)+iy+2;
            cS = ix*(Ny+1)+iy;
            cE = (ix+1)*(Ny+1) + iy+1;
            cW = (ix-1)*(Ny+1) + iy+1;
            

            
            row0 = [r,r,r,r,r];
            col0 = [r,cN,cE,cS,cW];
            
            for ii=1:5
                if(abs(star(ii)) > 0)
                    row(count+1,1) = row0(ii); 
                    col(count+1,1) = col0(ii); 
                    data(count+1,1) = star(ii); 
                    count = count+1; 
                end
            end
            
            
            b(r,:) = rhs';
        else
            row(count+1,1) = ix*(Ny+1)+iy+1;
            col(count+1,1) = ix*(Ny+1)+iy+1;
            data(count+1,1) = 1;
            
            count = count + 1; 
        end
        
        
    end
end
row = row(1:count,1); 
col = col(1:count,1); 
data = data(1:count,1); 

N = (Nx+1)*(Ny+1); 
S = sparse(row,col,data,N,N); 
end    

function I_out = get_image_at_grid(ix,iy,I0,A_info,phase_index,color_flag)

% ix = 0,...,Nx
% iy = 0,...,Ny

% Iout: scalar image entry or 3x1 vector image entry in case of colored
% images

Nx = size(I0,2); 
Ny = size(I0,1); 


n=0;

if(color_flag == 0)
    I = 0; 
    for k=0:1
        for l=0:1
            jx = ix+k;
            jy = iy+l;
            
            if(jx>=1 && jx<=Nx && jy>=1 && jy<=Ny)
                
                if(A_info(jy,jx)==phase_index)
                    I=I+ I0(jy,jx);
                    n=n+1;
                end
            end
        end
    end
else
    I = [0;0;0]; 
    for k=0:1
        for l=0:1
            jx = ix+k;
            jy = iy+l;
            
            if(jx>=1 && jx<=Nx && jy>=1 && jy<=Ny)
                
                if(A_info(jy,jx)==phase_index)
                    for alpha = 1:3
                        I(alpha)=I(alpha)+ I0(jy,jx,alpha);
                    end
                    n=n+1;
                end
            end
        end
    end
end


if(n>0)
    I=I/n;
end
I_out = I; 
end