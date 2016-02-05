%-------------------------------------------------------------------------%
% Script:             main
%
% Description:        Image segmentation with a parametric method using
%                     active contours (Chan Vese) with topology changes, 
%                     software for curves in the plane and curves on 
%                     surfaces
%
%-------------------------------------------------------------------------%

clear all
close all

feature accel on

%% Initial preparations, initial curve


% Set dimension
dimension = 2;   % 2: 2-dim., planar image, curve in R^2
                 % 3: Image and curves on a 2-dim. surface in R^3


% Read init file (which calls parameter file), 
% creation of the structures: config, Image, Surface, initial_information
init;

% Get image
fprintf('Preparations:\nGet image\n');
[Image,Surface] = get_image(config,Image,Surface);


% Plot draw image
figure(config.plot.f1)
draw_image(config,Image,Surface);



% Get color flag (dependent on size of image data)
color_flag = 0;
if(config.dimension==2)
    if(size(size(Image.data),2)==3)  % Nx x Ny x d  (d=1 gray, d=3 color)
        color_flag = 1;
    end
    
else
    if(size(Image.data,2)==3)   % Nsimp x d, (d=1 gray, d=3 color)
        color_flag = 1;
    end
    
end


if(config.tracking == 0)     % case: processing of a single image
    % Get initial curve(s)
    fprintf('Get initial curve and neighbor info\n');
    Gamma = get_initial_curve(config,initial_information,Surface);
    Gamma = get_initial_neighbor_info(Gamma);
    
    
    % Plot initial curve
    figure(config.plot.f1)
    draw_image(config,Image,Surface);
    plot_curves(Gamma,Image);
    
    
    % Calc length, perform global refinement / coarsening if necessary
    fprintf('Get length info\n');
    Gamma = get_length_info(Gamma,config,Image,Surface);

    
    % Update X_old (important if refinement or coarsening has been performed)
    Gamma.X_old = Gamma.X;
    
    
    % Get initial regions
    fprintf('Get initial regions\n');
    Gamma = calc_normal_field(Gamma,Surface);
    

    if(dimension == 2)
        % 2D image segmentation   
        init_center = [];
        init_radius = [];
        for i=1:size(initial_information.center,3)
            init_center = [init_center, initial_information.center(:,:,i)];
            init_radius = [init_radius, initial_information.radius(:,:,i)];
        end
        if(color_flag == 0)
            I0 = double(Image.data);
        else
            I0 = double([Image.data(:,:,1), Image.data(:,:,2), Image.data(:,:,3)]);
        end
        [Omega.A_info,Omega.I_info,Omega.n_info,Omega.coeffs] = get_initial_regions_2d(init_center,init_radius,Gamma.orient,[Image.sizes, size(initial_information.center,3), size(initial_information.center,2),config.method.color],I0);
    else
        % segmentation of images on 2-dim surfaces
         Omega  = get_initial_region_info_geodesic(Gamma,initial_information,Surface,Image,config);
    end
        
    
else
    % Load initial data for tracking
    [Gamma,Omega] = get_initial_data_for_tracking(config,Image);

end


% Draw image
draw_image(config,Image,Surface);
% Plot initial curve
plot_curves(Gamma,Image);




% Adapt lambda if wrong dimensions
config = adapt_lambda(config,Omega,color_flag);
% Adapt config.length.tol
if(size(config.length.tol,2)==1) 
    config.length.tol = [config.length.tol, config.length.tol]; 
end



%% Main loop
Prev = zeros(0,config.dimension);

k = 1;      % Main step counter
k_sub = 0;  % Sub step counter

repeat = 0; % Repeat flag (a step can be repeated if a topology change is detected, in this case the main step is divided in substeps)

ktrack = 1;
N_images = max([1,config.tracking]);


node_count = zeros(N_images,1); 

while(ktrack <=N_images)
    config.trackingnr = ktrack;
    if(config.tracking > 0)
        fprintf('Get image\n');
        [Image,Surface] = get_image(config,Image,Surface);
    end
    
    
    fprintf('\n---------------------------------------------\nImage %d/%d:\n', ktrack, N_images);
    
    
    while(k<=config.K)
        if(mod(k,config.plot.nr1)==0)
            fprintf('\n---------------------------------------------\n %d\n', k);
        end
        
        % Set time step size
        delta_t = config.timestep;

        if(k-config.top_check.k_top_last< config.top_check.kwait)
            repeat = 1;
            if(k_sub == 0)
                k_sub = 1;
            end
        end

        
        if(repeat>0)
            delta_t = delta_t/config.top_check.n_sub;
        end
        
        % Set small time step for last 10 k's
        if(config.K-k<10 && config.K>30)
            delta_t = delta_t/10;
        end
       
        % Update Gamma.X_old
        Gamma.X_old = Gamma.X;
        
        % Calculate normal vector field Gamma.nu
        Gamma = calc_normal_field(Gamma,Surface);
        
        % Get region info
        Omega = get_region_info(Gamma,Omega,config,Image,Surface);
        
        
        % Compute internal and external energy for automatic setting of sigma - implemented and used only for planar images! 
        if(dimension == 2 && config.method.sigma_compute > 0 && (k==1 || mod(k,config.method.sigma_compute)==0) && (repeat == 0 || k_sub==1))
            [Eint,Eext] = compute_energies(config,Gamma,Image,Omega,Surface);
            config = set_sigma(config,Eint,Eext,config.method.sigma_factor); 
            Omega.energies.int = Eint;
            Omega.energies.ext = Eext; 
        end

        
        % Calculate right hand side of the linear system
        b = calc_right_hand_side_new(Gamma,Omega,Surface,config,Image); % umfpack
        
        
        
        if(dimension == 2)
            [Gamma,S,P,M,b2] = solve_umfpack(Gamma,b,delta_t,config,Image);
        else
            Gamma = solve_umfpack_geodesic(Gamma,b,delta_t,config,repeat);
        end
      
        
        % Update Gamma.X
        Gamma.X = Gamma.X_old + Gamma.delta_X;
        if(config.dimension==2)
            Gamma  = image_limit_check(Gamma,Image);
        end
        
        % Look for topological changes
        check0 = 0;
        if(config.top_check.flag == 1)
            % Detect topological change
            [split,merge,triple,boundary,info0,mark1] = detect_top_change(Gamma,Image,Surface,Prev,config,repeat);

            N_topchange = size(split,1)+size(merge,1)+size(triple,1)+size(boundary,1);
            
            % Set repeat flag
            if(N_topchange > 0 && repeat == 0)
                fprintf('Detected topological change, repeat step with smaller step size\n');
                repeat = 0.5;
                k_sub = 0;
            else
                if(N_topchange>0)
                    fprintf('Top change detected, perform top change now\n');
                    repeat = 2;
                end
            end
            % Set Prev to 0 if no top changes has happenend and repeat == 0
            if(repeat == 0 && N_topchange == 0)
                Prev = zeros(0,config.dimension);
            end
            
            
        end
        
        % Perform topological change
        if(repeat == 2)
            % Perform top change
            fprintf('Performing topological change...\n');
            [Gamma,Prev] = perform_top_change(split,merge,triple,boundary,Gamma,Image,Surface,config,Prev);       
            fprintf('Finished topological change.\n');
            repeat = 3;
            
            config.top_check.k_top_last = k;
            
        end
        
        
        % Calculate length of the curves
         if(repeat ~=0.5)  % Do not check this while you start repeating the step
            % Calculate length info and mark curves for deleting (if their length
            % is too small)
            
            Gamma = get_length_info(Gamma,config,Image,Surface);
            
            
            % Delete curves and corresponding nodes, reset J, X, neigh,
            % index_info and nr_curves
            if(sum(sum(Gamma.mark))>0)
                
                fprintf('Delete curve\n');
                % Delete too small curves
                Gamma = delete_curves(Gamma);
                Gamma = calc_normal_field(Gamma,Surface);
            end
         end
        
        
        % Update closest simplex for a node point on a triangulated surface
        % (dim = 3 only, i.e. curves on surfaces)
        if(config.dimension==3)
            if(Image.flag == 7)
                project_to_surface_torus
                Gamma = update_closest_simp(Gamma,Surface);
                project_to_surface_torus
            else
                if(Image.flag ==5 || Image.flag == 6)
                    project_to_surface_sphere;
                    Gamma = update_closest_simp(Gamma,Surface);
                    project_to_surface_sphere;
                else
                    Gamma = update_closest_simp(Gamma,Surface);
                end
            end
        end       
        
        
        % Print information and show/plot curves and image
        if(mod(k,config.plot.nr1)==0)
            % Plot curves
            draw_image(config,Image,Surface);
            figure(config.plot.f1);
            plot_curves(Gamma,Image);
            
            % Draw other viewing direction in case some surfaces 
            if(config.dimension == 3)
                if(Surface.flag >= 5)
                    figure(config.plot.f2);
                    plot_curves(Gamma,Image); 
                    figure(config.plot.f3); 
                    plot_curves(Gamma,Image); 
                end
            end
            
            % Print information
            fprintf('k=%d, nr of curves =',k);
            for ll=1:size(Gamma.index_info,1)
                fprintf(' %d  ', Gamma.nr_curves(ll,1));
            end
            fprintf('nr of nodes %d\n', size(Gamma.X,1)); 
            fprintf('\n'); 

            

            fprintf('delta_t = %2.3e, repeat = %d\n\n',delta_t, repeat);
            
        end
        
        % End sub-loop if k_sub == n_sub
        if(k_sub == config.top_check.n_sub)
            repeat = 0;
            k_sub  = 0;
        end
 
        
        % End step 
        if(repeat == 0)
            % Update k
            k=k+1;
            % Update age of curves
            for ll=1:size(Gamma.index_info,1)
                Gamma.age(ll,1:Gamma.nr_curves(ll,1)) = Gamma.age(ll,1:Gamma.nr_curves(ll,1))+1;
            end
            
        else
            if(repeat == 0.5)  % intermediate value for repeat
                repeat  = 1; 
            end
            
            if(k_sub == 0) % If main loop should be divided in sub-steps (e.g. close to a topology change)
                Gamma.X = Gamma.X_old;  % Reset X to X_old repeat the last step and divide it in 10 sub-steps
            end
            % Update k_sub
            k_sub = k_sub + 1;
            
        end
        
        
    end
    
    ktrack = ktrack + 1;
    k=1;
end



% After final step
draw_image(config,Image,Surface);
plot_curves(Gamma,Image);

fprintf('End of segmentation\nPress Enter to continue with post image denoising.\n'); 
pause



%% After final step perform image smoothing
if(config.dimension == 2)
   Image_new = image_diffusion_and_denoising_2d_sparse(Image,Omega,config);
else
    u = diffusion_eq_geodesic2(Surface,Image,Omega,1/config.param_diffusion);
    Image_new = Image; 
    Image_new.data = u; 
    Image_new.draw = u; 
    draw_image(config,Image_new,Surface); 
end




