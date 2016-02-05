function draw_cross_sections(Gamma,Omega,y_const,z_const)



Nx = size(Omega.data,2); 
Ny = size(Omega.data,1); 
Nz = size(Omega.data,3); 

cmap = zeros(255,3); 
for i=1:255
    cmap(i+1,:) = i/255;
end


%% z = const
z_const = Nz+1-z_const;

for i0=1:size(z_const,1)
    
    z=z_const(i0);
    points = compute_cross_sections(Gamma,z,3);
    
    I0 = Omega.data(:,:,z);
    
    figure(4+i0)
    clf
    image(I0);
    axis image
    colormap(cmap);
    
    % Adjust y-coord of points
    points_plot = points;
    for i=1:size(points_plot,1)
        points_plot(i,2) = Ny+1-points_plot(i,2);
    end
    hold on
    plot(points_plot(:,1), points_plot(:,2), 'b.');
    pause(0.1);
end

%% y = const
for i0=1:size(y_const,1)
    
    y = y_const(i0);
    points = compute_cross_sections(Gamma,y,2);
    
    I1 = uint8(zeros(Nz,Nx));
    for j=1:Nz
        for jj=1:Nx
            I1(j,jj) = Omega.data(Ny+1-y,jj,Nz+1-j);
        end
    end
    
    figure(7+i0)
    clf
    image(I1);
    colormap(cmap);
    axis image
    
    % Adjust z-coord of points
    points_plot = [points(:,1), points(:,3)];
    for i=1:size(points_plot,1)
        points_plot(i,2) = Nz+1-points_plot(i,2); 
    end
    hold on
    plot(points_plot(:,1), points_plot(:,2), 'b.');
    pause(0.1);
end
end