function plot_curves(Gamma,Image)

dim = size(Gamma.X,2);

Nmain = size(Gamma.index_info,1);

for k=1:Nmain
    
    for l=1:Gamma.nr_curves(k,1)
       
        
        boundary_points = zeros(2,dim);
        boundary_flag   = 0;
        
        X_plot=zeros(size(Gamma.X,1),dim);
        X_plot(1,:) = Gamma.X(Gamma.index_info(k,l),:);
        
        if(Gamma.neigh(Gamma.index_info(k,l),1)<0)
            boundary_points(1,:)=Gamma.X(Gamma.index_info(k,l),:);
            boundary_flag = 1;
        end
        
        i=1;
        j=Gamma.index_info(k,l);
        while(j~=Gamma.index_info(k,l) || i==1)
            if(Gamma.neigh(j,2)>0)
                X_plot(i+1,:)=Gamma.X(Gamma.neigh(j,2),:);
                i = i+1;
                j = Gamma.neigh(j,2);
            else
                boundary_points(2,:)=Gamma.X(j,:);
                j = Gamma.index_info(k,l); %stop loop

               
            end
            
        end
        if(dim ==2)
            X_plot(:,2)=Image.sizes(2)+1-X_plot(:,2);
        end
        
        

        
        switch dim
            case 2
                % Plot dots
                mp=plot(X_plot(1:i,1),X_plot(1:i,2), '.');
                hold on
            case 3
                str = 'b.';
                mp = plot3(X_plot(1:i,1),X_plot(1:i,2),X_plot(1:i,3),str);
                hold on
        end
        
        switch mod(k,15)
            case 1
                set(mp, 'Color', [1, 0, 0]);
            case 2
                set(mp, 'Color', [0, 0, 1]);
            case 3
                set(mp, 'Color', [0, 1, 0]);
            case 4
                set(mp, 'Color', [0, 1, 1]);
            case 5
                set(mp, 'Color', [1, 0, 1]);
            case 6
                set(mp, 'Color', [1, 1, 0]);
            case 7
                set(mp, 'Color', [1, 0.5, 0.5]);
            case 8
                set(mp, 'Color', [0.5, 1, 0.5]);
            case 9
                set(mp, 'Color', [0.5, 0.5, 1]);
            case 10
                set(mp, 'Color', [1, 0.5, 0]);
            case 11
                set(mp, 'Color', [1, 0, 0.5]);
            case 12
                set(mp, 'Color', [0.5, 1, 0]);
            case 13
                set(mp, 'Color', [0, 1, 0.5]);
            case 14
                set(mp, 'Color', [0, 0.5, 1]);
            otherwise
                set(mp, 'Color', [0.5, 0, 1]);
        end
        
    end
end
pause(0.02)
end