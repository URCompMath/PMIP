%-------------------------------------------------------------------------%
% Scipt:              main
%
% Description:        Implementation of parametric method for 3D image 
%                     segmentation
%                     (I)   V_n = X_t * nu = sigma kappa + F,
%                     (II)  kappa * nu = X_ss
%                     
%                     F is an external forcing term accoring to the
%                     piecewise constant Mumford-Shah problem
%
% Author:             Heike Benninghoff
% Date:               2012-2015
%-------------------------------------------------------------------------%

clear all
close all

feature accel on
profile on 


% Load parameters
init_parameters

% Load initial data
Gamma = load_surface(Gamma,config);

% Plot initial surface
addpos = [1865,0,0,0];
figure('Position', [ -1895   564    592    420] + addpos);
figure('Position', [ - 545   564    592    420] + addpos);
plot_surface(Gamma.X,Gamma.simplices,config.image.flag);
figure('Position', [ -1895    43    592    420] + addpos);
figure('Position', [ - 545    43    592    420] + addpos);
pause(0.1);




% Intialize regions (Omega)
if(config.method.pde_flag == 1)
    % Initial region information for image segmentation
    fprintf('Get initial region info\n');
    Omega = get_initial_region_info(Gamma,config,initial_info);
    for k=1:Omega.nr_regions
        fprintf('Region %d:   I_info = %3.2f  n_info=%d  coeffs = %3.6f\n', k, Omega.I_info(k,:), Omega.n_info(k,1), Omega.coeffs(k,:));
    end
else
    % E.g. mean curvature flow
    Omega = [];
end



% Main part
fprintf('Main loop\n');
k=1;
config.topchange.counter = config.topchange.countermin;

Gamma.time = 0;
Gamma.delta_t = config.timestepsize.delta_t_init;
config.T = config.K * Gamma.delta_t;
    

if(config.image.flag >= 6) % Medical CT Data, Additional Figure for drawing 2D cross-sections    
    figure('Position', [  0   400    500    360]);
    figure('Position', [450   400    500    360]);
    figure('Position', [900   400    500    360]);
    figure('Position', [  0    10    500    360]);
    figure('Position', [450    10    500    360]);
    figure('Position', [900    10    500    360]);
   
    draw_cross_sections(Gamma,Omega,config.cross.y_const,config.cross.z_const); 
end

sigma_t0 = max(config.mesh.sigma_t);
compute_external = 1; 
count_external = 0;

while(Gamma.time <= config.T)
    
    % Increase lambda for last time steps
    if(k> 500 && config.image.flag == 6) % Change lambda for lung example
        config.method.lambda = 1000;
        config.timestepsize.tol_DXn_max = 10; 
    else
        if(k> 500 && config.image.flag == 7)
            config.method.lambda = 50;
        end
    end
    
    % Set time step
    Gamma.delta_t = config.timestepsize.delta_t_init / config.timestepsize.factor_delta_t;
    
    % Set time
    Gamma.time = Gamma.time + Gamma.delta_t; 
   
    fprintf('\n\n-------------------------------------------------------------------\nTime Step %d   delta_t = %2.8f\nTime = %2.8f  Time end = %2.8f\n\n', k, Gamma.delta_t, Gamma.time, config.T); 
    
    
    % Improve mesh
    if(config.mesh.improve_flag && mod(k,10)==0)
        fprintf('Improve mesh\n'); 
        max_nr_simp = 9;
        max_locations = 30;
        Gamma = improve_mesh2_general(Gamma,config, max_nr_simp, max_locations);
    end
    
    % Compute normal and area
    [Gamma,mark_simp,mark_delete] = compute_normal_and_area(Gamma,config); 
    
    
    % Refine or delete marked simplices
    if(config.mesh.refinecoarsen.flag)
        mesh_flag = 0; 
        [X_new,simplices_new,edges_new,involved_simplices] = refine_marked_simplices(Gamma.X,Gamma.simplices,Gamma.edges,mark_simp);
        if(sum(mark_simp)>0)
             Gamma.X = X_new; 
             Gamma.simplices  = simplices_new;
             Gamma.edges = edges_new; 
             [Gamma, mark_delete] = compute_normal_and_area_delete_only(Gamma,config);
        end
        
        
      if(config.mesh.delete_flag)
        [X_new,simplices_new,edges_new] = delete_too_small_simplices(Gamma.X,Gamma.simplices,Gamma.edges,mark_delete,involved_simplices); 
        if(sum(mark_delete)>0)
            Gamma = area_compute_and_zero_check(Gamma,X_new,simplices_new,edges_new);
        end
      end
        
    end
   
    % Compute regions
    if(config.method.pde_flag == 1)
        fprintf('Get region info and compute coeffs\n');
        Omega = get_region_info(Gamma,Omega,config);
        
        for l=1:Omega.nr_regions
            fprintf('Region %d:   I_info = %3.2f  n_info=%d  coeffs = %3.6f\n',l, Omega.I_info(l,:), Omega.n_info(l,1), Omega.coeffs(l,:)); 
        end
        fprintf('\n'); 
    end
    
    
    if(config.mesh.flag == 2 && config.topchange.counter >= config.topchange.countermin)
            config.mesh.sigma_t = get_sigma(Gamma,sigma_t0,config);
    else
        % Default (color for plot surface)
        config.mesh.sigma_t = sigma_t0 * ones(size(Gamma.X,1),1);
    end

    
    if(compute_external)
        fprintf('compute_external = 1, solve Vn = sigma kappa + F\n'); 
    else
        fprintf('compute_external = 0, solve Vn = sigma kappa\ncount_external = %d\n', count_external); 
        count_external = count_external + 1;
        if(count_external >= 0)
            compute_external = 1;
        end
    end
    
    % Solve linear system
    fprintf('Solve linear system\n'); 
    Gamma = solve_umfpack(config,Gamma,Omega,compute_external);
    

    % Time step control
    repeat = 0; 
    if(config.timestepsize.flag)
        fprintf('\nTime step size control:\n'); 
        Gamma = compute_omega(Gamma);
        DXn = compute_deltaX_normal(Gamma);
        DXnmax = max(abs(DXn));
        fprintf('max abs dXn  = %2.5f   tol_dXn_min = %2.5f, tol_dXn_max = %2.5f\n', DXnmax,config.timestepsize.tol_DXn_min,config.timestepsize.tol_DXn_max);
        
        if(DXnmax > config.timestepsize.tol_DXn_max)
            config.timestepsize.factor_delta_t = config.timestepsize.factor_delta_t * config.timestepsize.decrease_factor;
            repeat = 1;
            Gamma.time = Gamma.time - Gamma.delta_t;
            fprintf('Repeat and decrease time step size.\n');
        else
            if(DXnmax < config.timestepsize.tol_DXn_min)
                config.timestepsize.factor_delta_t = config.timestepsize.factor_delta_t / config.timestepsize.decrease_factor;
                repeat = 1;
                Gamma.time = Gamma.time - Gamma.delta_t;
                fprintf('Repeat and increase time step size.\n');
            end
        end
    end

    
    % Proceed if repeat == 0
    if(repeat == 0)
        
        % Update X (delta_X = X - X_old)
        Gamma.X = Gamma.delta_X + Gamma.X;
        
        % Plot surface before mesh changes
        if(size(Gamma.X,1)>0)
            figure(1)
            plot_surface(Gamma.X,Gamma.simplices,config.image.flag);
        end
                
        
        
        % Topological changes
        if(config.topchange.flag==1)

            fprintf('\nCheck for top changes\n');
            [Gamma, top_change] = topological_changes_four_cases(Gamma,config);
            
            if(top_change > 0)
                config.mesh.sigma_t = zeros(size(Gamma.X,1),1); 
                config.topchange.counter = 0;
            end


        end

        
        % Calc area of single surfaces and mark for deletion
        [Gamma.area_info,mark] = get_area_info(Gamma,config);
        
        % Delete surfaces if it is marked for deletion (if area < area_tol)
        Gamma = delete_surfaces(Gamma,mark,config); 

        
        % Plot surface
        if(size(Gamma.X,1)>0)
            figure(2)
            plot_surface2(Gamma.X,Gamma.simplices,config);
            figure(3)
            plot_surface_cross(Gamma.X,Gamma.simplices,0,0,config.image.flag)
            figure(4)
            plot_surface_with_viewing_angle(Gamma.X,Gamma.simplices,0,90,config.image.flag)
            
            if(config.image.flag >= 6)
                draw_cross_sections(Gamma,Omega,config.cross.y_const,config.cross.z_const);
            end
        end
        
        
        % Write data
        fprintf('\nResults:\n'); 
        for l=1:Gamma.nr_surfaces
            fprintf('k=%d, J=%d  surface nr = %d , area = %2.5f\n', k, size(Gamma.X,1),l, Gamma.area_info(l));
            
        end
               
        % Update
        k=k+1;
        config.topchange.counter = config.topchange.counter + 1; 

        if(Gamma.nr_surfaces==0)
            k=config.K+1;
            Gamma.time = config.T+1; 
        end

        
        
    end
    
end

profile viewer