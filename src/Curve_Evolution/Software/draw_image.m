function draw_image(config,Image,Surface)

figure(config.plot.f1);
clf;
switch config.dimension
    case 2   % 2: 2-dim., planar image, curve in R^2
        if(size(size(Image.data),2)==2)  % scalar image, need colormap
            image(Image.draw);
            colormap(Image.cmap);
            axis image
        else
            image(Image.draw);           % vector valued image, color image
            axis image
        end
        hold on
        
        
    case 3   % 3: Image and curves on a 2-dim. surface in R^3
        % Get color data
        if(size(Image.data,2)==1)
            cdata = Image.draw;  % gray image
        else
            if(max(max(Image.draw))> 1)
                cdata = Image.draw/255; % scale to [0,1];
            else
                cdata = Image.draw;
            end
        end
        
        p=patch('Vertices',Surface.nodes,'Faces',Surface.tri);
        set(p,'FaceColor','flat','FaceVertexCData',cdata,'edgecolor','none','CDataMapping','direct');
        
        set(p,'FaceLighting','gouraud',...
            'AmbientStrength',.5,'DiffuseStrength',.5,...
            'SpecularStrength',.2,'SpecularExponent',25,...
            'SpecularColorReflectance', 1);
        set(p,'FaceAlpha',1 );
        
        if(size(Image.data,2)==1) % colormap for gray images
            colormap(Image.cmap);
        end
        
        axis equal
        
        % Change viewing angle
        switch Surface.flag

            case 1
                view(0,90)
            case 5
                view(75,15)
            case 6
                view(75,15)
            case 7
                view(75,40); 
        end
        xlabel('x');
        ylabel('y');
        zlabel('z');
        hold on
        
        
        
        if(Surface.flag >= 5) % Earth or Torus, plot 2 additional views
            figure(config.plot.f3);
            clf;
            p=patch('Vertices',Surface.nodes,'Faces',Surface.tri);
            set(p,'FaceColor','flat','FaceVertexCData',cdata,'edgecolor','none','CDataMapping','direct');
            
            set(p,'FaceLighting','gouraud',...
                'AmbientStrength',.5,'DiffuseStrength',.5,...
                'SpecularStrength',.2,'SpecularExponent',25,...
                'SpecularColorReflectance', 1);
            set(p,'FaceAlpha',1 );
            
            if(size(Image.data,2)==1) % colormap for gray images
                colormap(Image.cmap);
            end
            
            lighting flat
            axis equal
            
            view(195,15);
            
            if(Surface.flag == 7)
                view(195,40); 
            end
            xlabel('x');
            ylabel('y');
            zlabel('z');
            hold on
            
            
            
            figure(config.plot.f4);
            clf;
            p=patch('Vertices',Surface.nodes,'Faces',Surface.tri);
            set(p,'FaceColor','flat','FaceVertexCData',cdata,'edgecolor','none','CDataMapping','direct');
            
            set(p,'FaceLighting','gouraud',...
                'AmbientStrength',.5,'DiffuseStrength',.5,...
                'SpecularStrength',.2,'SpecularExponent',25,...
                'SpecularColorReflectance', 1);
            set(p,'FaceAlpha',1 );
            
            if(size(Image.data,2)==1) % colormap for gray images
                colormap(Image.cmap);
            end
            
            lighting flat
            axis equal
            
            view(315,15);
            if(Surface.flag == 6)
                view(315,-25);
            else
                if(Surface.flag == 7)
                    view(315,40); 
                    
                end
            end
            
            xlabel('x');
            ylabel('y');
            zlabel('z');
            hold on
        end
        
end






pause(0.01);




end