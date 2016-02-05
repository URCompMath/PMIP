function plot_surface2(X,simplices,config)

% Plot surface
N = size(simplices,1);

tri = zeros(N,3);
for i=1:N
    tri(i,:)=simplices{i,1}.nodes;
end


C = zeros(size(X));
for i=1:size(X,1)
    if(config.mesh.sigma_t(i)==0)
        C(i,:) = [1 0 0];
    else
        C(i,:) = [0 1 0];
    end
end



clf

p = patch('Vertices', X, 'Faces', tri);
set(p,'FaceColor','flat','FaceVertexCData',C,'edgecolor', [0 0 0] ,'CDataMapping','direct');

set(p,'FaceLighting','gouraud',...
    'AmbientStrength',.5,'DiffuseStrength',.5,...
    'SpecularStrength',.2,'SpecularExponent',25,...
    'SpecularColorReflectance', 1);
set(p,'FaceAlpha',1 );

% Default figure axis sizes
axis([-2.5 2.5 -1.7 1.7 -1.7 1.7])

image_flag = config.image.flag;        
if(image_flag == 3) % Merging
    axis([-2.5 2.5 -1.7 1.7 -1.7 1.7]/2)
else
    if(image_flag == 6) % Lung Segm. TCIA
        axis([-1 446 -1 310 -1 276]);
    else
        if(image_flag == 7) % CT data UKR Lung+Heart Segm.
            axis([-1 129 -1 129 -1 142]);
        else
            if(image_flag == 8) % CT data UKR Abdominal Region Segm.
                axis([-1 129 -1 129 -1 81]);
            else
                if(image_flag == 9)  % CT data UKR Lung Segm., Splitting
                    axis([-1 129 -1 129 -1 281]);
                end
            end
        end
    end
end

view(-5,45); 




grid on
pause(0.1)
end