function [config,Image,Surface,initial_information] = init_parameters_geodesic()
% PARAMETER FILE

% Image flag
image_flag = 1;   
                    % 1: bunny with black circles
                    % 2: face001
                    % 3: face052
                    % 4: face053
                    % 5: Earth - longwave radiation
                    % 6: Earth - net radiation
                    % 7: Torus Objects

     
% Surface flag (here identical to image flag)
surface_flag = image_flag;  
                    % 1: bunny surface
                    % 2: face001 
                    % 3: face052 
                    % 4: face053
                    % 5: Earth
                    % 6: Earth 
                    % 7: Torus Objects
                    
%% Set the parameters

%% Default parameters
% (will be partly overwritten dependent on image_flag and surface_flag)
                    
% Analytical surface or triangulated surface (incl. image data)?
% image segmentation examples: here only triangulated surfaces, parts of
% the software can also be used for analytical surfaces (e.g. mean
% curvature flow of curves on a sphere. Here however image data is always
% given on triangulated images
triangulation = 1;  % 0: analytical surface, 1: triangulated surface


                    
% Image Segmentation Method 
im_seg_method = 1;   % 1: Chan Vese with parametric method incl. multiple
% phases, incl. triple junctions and topology changes
% (parameter im_seg_method can be used to implement other segmentation
% methods like Geodesic Active Contours here, currently only Chan-Vese)


% Lmax and Lmin, parameters for global refinement/ coarsening
Lmax =  1; 
Lmin =  0.2;  
Nmin = 10;  

% Use different length tolerance for new curves (new small curves should
% not be immediately deleted, time for growing --> age_min time steps)
age_min = 5; 

% Width of narrow band where the regions are updated around Gamma
width_band = 5; 

% Substeps if main step is divided
n_sub = 10; 
k_top_last = -100; 
kwait = 2;

% Grid size top change
amin = 1; 
amax = 3; 

% Initialization of parameters 
sigma = 1; 
lambda = 1; 
lambda_partial = 0.8*sigma;
mu = 1;

% Default color method
color_method = 0;           
% 0: gray-scaled, no color 
% 1: RGB
% 2: CB
% 3: HSV

% Default tracking flag
tracking = 0;

% Tolerance used in solving linear system with cg-method
tol_std = 1e-4;   
tol_extra = 1e-4;

% Band width for Prev in detect_top_change
a_Prev = 0.5; 

% Setting of sigma
sigma_factor = 0.2;
sigma_compute = 0;
sigma_min = 30;
sigma_step = 5;

% Surface tol (for is_inside decisions)
surface_tol = 0.005;
level_max = 10;

%% Image / Surface dependent parameters

switch image_flag

 
    case 1    % bunny, 3 circles
        % Image sizes
        size_x     = 158;
        size_y     = 157;
        size_z     = 122;
        % Number of nodes (initial number)
        J = 30; 
        % Number of timesteps
        K = 160; 
        % Main time step size
        Delta_t = 0.01;
        % Parameter sigma (weight for curvature)
        sigma = 2;
        % Parameter lambda (weight for external force)
        lambda = 50; 
        % Parameter diffusion  - mu Laplace u + u = u_0
        param_diffusion = 1; 
        % Plot number 
        plot_nr = 10; 
        plot_nr2 = 10; 
        % Check topological changes? 0: no, 1:yes
        top_check = 1;
        % Length tolerance
        length_tol = 2*pi*0.02; 
        % Coarsening / Refinement parameter
        Lmax =  0.25; 
        Lmin =  0.025;
        
    case 2    % face001 with color
        % Image sizes
        size_x     = 176;
        size_y     = 209;
        size_z     = 143;
        % Number of nodes (initial number)
        J = 150;
        % Number of timesteps
        K = 70; 
        % Main time step size
        Delta_t = 0.01;
        % Parameter sigma (weight for curvature)
        sigma = 0.25; 
        % Parameter lambda (weight for external force)
        lambda = [200,20]; 
        % Color method
        color_method = 2;  % 2: CB
        % Parameter diffusion  - mu Laplace u + u = u_0
        param_diffusion = 1; 
        % Plot number 
        plot_nr = 1; 
        plot_nr2 = 1; 
        % Check topological changes? 0: no, 1:yes
        top_check = 0;
        % Length tolerance
        length_tol = 2*pi*0.02; 
        % Coarsening / Refinement parameter
        Lmax =  0.15; 
        Lmin =  0.025; 
        % Width band
        width_band = 4;    
        % Max search levels for closest simplex
        level_max = 20;
        
   case 3    % face052 with color
        % Image sizes
        size_x     = 184;
        size_y     = 192;
        size_z     = 146;
        % Number of nodes (initial number)
        J = 80;
        % Number of timesteps
        K = 50; 
        % Main time step size
        Delta_t = 0.002;
        % Parameter sigma (weight for curvature)
        sigma = 0.25; 
        % Parameter lambda (weight for external force)
        lambda = [200,20]; 
        % Color method
        color_method = 2;  % 2: CB
        % Parameter diffusion  - mu Laplace u + u = u_0
        param_diffusion = 1; 
        % Plot number 
        plot_nr = 1; 
        plot_nr2 = 1; 
        % Check topological changes? 0: no, 1:yes
        top_check = 0;
        % Length tolerance
        length_tol = 2*pi*0.02; 
        % Coarsening / Refinement parameter
        Lmax =  0.1; 
        Lmin =  0.01;      
        % Width band
        width_band = 4;
        % Max search levels for closest simplex
        level_max = 20;
        
    case 4    % face053 with color
        % Image sizes
        size_x     = 191;
        size_y     = 220;
        size_z     = 155;
        % Number of nodes (initial number)
        J = 80;
        % Number of timesteps
        K = 60; 
        % Main time step size
        Delta_t = 0.002;
        % Parameter sigma (weight for curvature)
        sigma = 0.25; 
        % Parameter lambda (weight for external force)
        lambda = [200,20]; 
        % Color method
        color_method = 2;  % 2: CB
        % Parameter diffusion  - mu Laplace u + u = u_0
        param_diffusion = 1; 
        % Plot number 
        plot_nr = 1; 
        plot_nr2 = 1; 
        % Check topological changes? 0: no, 1:yes
        top_check = 1;
        % Length tolerance
        length_tol = 2*pi*0.02; 
        % Coarsening / Refinement parameter
        Lmax =  0.15; 
        Lmin =  0.01;      
        % Width band
        width_band = 4;
        % Max search levels for closest simplex
        level_max = 20;

    case 5  % Earth - longwave radiation
        % Image sizes
        size_x     = -1;
        size_y     = -1;
        size_z     = -1;
        % Number of nodes (initial number)
        J = 6*20;
        % Number of timesteps
        K = 80; 
        % Main time step size
        Delta_t = 0.001;
        % Parameter sigma (weight for curvature)
        sigma = 1; 
        % Parameter lambda (weight for external force)
        lambda = 50; 
        % Color method
        color_method = 1; % RGB
        % Parameter diffusion  - mu Laplace u + u = u_0
        param_diffusion = 1; 
        % Plot number 
        plot_nr = 1; 
        plot_nr2 = 1; 
        % Check topological changes? 0: no, 1:yes
        top_check = 1; 
        % Length tolerance
        length_tol = 2*pi*0.002; 
        % Coarsening / Refinement parameter
        Lmax =  0.08; 
        Lmin =  0.005;      
        % Width band
        width_band = 2;  
        % Max search levels for closest simplex
        level_max = 30;

    case 6  % Earth - net radiation
        % Image sizes
        size_x     = -1;
        size_y     = -1;
        size_z     = -1;
        % Number of nodes (initial number)
        J = 120;
        % Number of timesteps
        K = 200; 
        % Main time step size
        Delta_t = 0.002; 
        % Parameter sigma (weight for curvature)
        sigma = 1;  
        % Parameter lambda (weight for external force)
        lambda = 300;
        % Color method
        color_method = 1; % RGB
        % Parameter diffusion  - mu Laplace u + u = u_0
        param_diffusion = 1; 
        % Plot number 
        plot_nr = 1; 
        plot_nr2 = 1; 
        % Check topological changes? 0: no, 1:yes
        top_check = 1; 
        % Length tolerance
        length_tol = 2*pi*0.002; 
        % Coarsening / Refinement parameter
        Lmax =  0.08; 
        Lmin =  0.005;      
        % Width band
        width_band = 2;  
        % Max search levels for closest simplex
        level_max = 30;
       
    case 7  % Torus Objects
        % Image sizes
        size_x     = -1;
        size_y     = -1;
        size_z     = -1;
        % Number of nodes (initial number)
        J =6*44; 
        % Number of timesteps
        K = 470; 
        % Main time step size
        Delta_t = 0.0001;
        % Parameter sigma (weight for curvature)
        sigma = 1; % 1; 
        % Parameter lambda (weight for external force)
        lambda = 20; 
        % Change color method
        color_method = 1; % 1: RGB
        % Parameter diffusion  - mu Laplace u + u = u_0
        param_diffusion = 1; 
        % Plot number 
        plot_nr = 1; 
        plot_nr2 = 1; 
        % Check topological changes? 0: no, 1:yes
        top_check = 1; 
        % Length tolerance
        length_tol = 2*pi*0.002; 
        % Coarsening / Refinement parameter
        Lmax =  0.06; 
        Lmin =  0.02;      
        % Width band
        width_band = 4; 
        % Max search levels for closest simplex
        level_max = 40;
        
end


switch surface_flag

    case 1 % bunny, 3 circles
        R = 1;             % not used
        r0 = 1;            % not used
        r = zeros(3,1,1);  % radius of initial curves
        r(:,1,1) = 1.8*[1;1;0]; 
        c = zeros(3,1,1);  % center of initial curves
        c(:,1,1) = [0;8;5.8]; 
        % Change default amin, amax
        amin = 0.05;
        amax = 0.3; 

    case 2 % face001
        R = 1;             % not used
        r0 = 1;            % not used
        r = zeros(3,1,1);  % radius of initial curve ---> mouth
        r(:,1,1) = 2*[2;0.8;0]; 
        c = zeros(3,1,1);  % center of initial curve
        c(:,1,1) = [-0.25;-3.2;10.8]; 
        % Change default amin, amax
        amin = 0.05;
        amax = 0.1; 
        
    case 3 % face052
        R = 1;             % not used
        r0 = 1;            % not used
        r = zeros(3,1,1);  % radius of initial curve
        r(:,1,1) = 1.5*[2;0.8;0]; 
        c = zeros(3,1,1);  % center of initial curve
        c(:,1,1) = [0;-3.1;13];        
        % Change default amin, amax
        amin = 0.05;
        amax = 0.3; 
        
    case 4 % face053
        R = 1;             % not used
        r0 = 1;            % not used
        r = zeros(3,1,1);  % radius of initial curve
        r(:,1,1) = 1.5*[2;1;0]; 
        c = zeros(3,1,1);  % center of initial curve
        c(:,1,1) = [0;-3.5;13.6];        
        % Change default amin, amax
        amin = 0.05;
        amax = 0.3; 
        
    case 5   % Earth - longwave radiation
        R = 1;             % not used
        r0 = 1;            % not used
        r = zeros(3,4,1);  % radius of initial curve 
        r(:,1,1) = [0;0.3;0.3]; 
        r(:,2,1) = [0.3;0;0.3];
        r(:,3,1) = [0;0.3;0.3]; 
        r(:,4,1) = [0.3;0;0.3]; 
        c = zeros(3,4,1);  % center of initial curve
        c(:,1,1) = [1;0;0];  
        c(:,2,1) = [0;1;0];
        c(:,3,1) = [-1;0;0];
        c(:,4,1) = [0;-1;0]; 
        % Change default amin, amax
        amin = 0.005;
        amax = 0.08;
        
    case 6   % Earth - net radiation
        R = 1;             % not used
        r0 = 1;            % not used
        r = 0.32* ones(3,2,2);  % radius of initial curve 
        c = zeros(3,2,2);  % center of initial curve
        c(:,1,1) = [0;0;1];
        c(:,2,1) = [0;0;-1];
        offset = 32*pi/180;
        beta = 0; 
        alpha = 0 + offset;
        c(:,1,2) = [cos(alpha)*cos(beta);sin(alpha)*cos(beta);sin(beta)];
        alpha = pi + offset;
        c(:,2,2) = [cos(alpha)*cos(beta);sin(alpha)*cos(beta);sin(beta)];
        % Change default amin, amax
        amin = 0.005;
        amax = 0.08;
      
    case 7   % Torus objects
        R = 1;             % not used
        r0 = 1;            % not used
        r = 0.4* ones(3,2,3);  % radius of initial curve 
        c = zeros(3,2,3);  % center of initial curve
        alpha = pi/6; 
        c(:,1,1) = 1.2* [cos(alpha); sin(alpha); 0] + 0.05*[0;0;1]; 
        alpha = alpha + pi/3; 
        c(:,1,2) = 1.2* [cos(alpha); sin(alpha); 0] + 0.05*[0;0;1]; 
        alpha = alpha + pi/3; 
        c(:,1,3) = 1.2* [cos(alpha); sin(alpha); 0] + 0.05*[0;0;1]; 
        alpha = alpha + pi/3; 
        c(:,2,1) = 1.2* [cos(alpha); sin(alpha); 0] + 0.05*[0;0;1]; 
        alpha = alpha + pi/3; 
        c(:,2,2) = 1.2* [cos(alpha); sin(alpha); 0] + 0.05*[0;0;1]; 
        alpha = alpha + pi/3; 
        c(:,2,3) = 1.2* [cos(alpha); sin(alpha); 0] + 0.05*[0;0;1]; 
        % Change default amin, amax
        amin = 0.005;
        amax = 0.02;
end



f1=figure;
set(f1, 'Position', [136   550   560   420]);

if(surface_flag >= 5)
    f2=figure;
    set(f2, 'Position', [712   550   560   420]);
    f3=figure;
    set(f3, 'Position', [136    50   560   420]);
end


%% Write parameters to structure

% Configuration parameters
config.dimension = 3; 
config.K = K; 
config.timestep = Delta_t; 

config.method.flag = im_seg_method;
config.method.sigma_factor = sigma_factor;
config.method.sigma_compute = sigma_compute;
config.method.sigma_min = sigma_min;
config.method.sigma_step = sigma_step;
config.method.sigma = sigma; 
config.method.lambda = lambda; 
config.method.lambda_partial = lambda_partial;
config.method.mu = mu; 
config.method.color = color_method;
config.param_diffusion = param_diffusion; 

config.tol.std = tol_std;
config.tol.extra = tol_extra;

config.length.tol = length_tol; 
config.length.Lmax = Lmax;
config.length.Lmin = Lmin;
config.length.Nmin = Nmin; 
config.length.age_min = age_min; 

config.width_band = width_band; 


config.top_check.flag = top_check; 
config.top_check.n_sub = n_sub;
config.top_check.k_top_last = k_top_last;
config.top_check.kwait = kwait; 
config.top_check.amin = amin;
config.top_check.amax = amax;
config.top_check.a_Prev = a_Prev;

config.plot.nr1 = plot_nr;
config.plot.nr2 = plot_nr2; 
config.plot.f1  = f1; 
if(surface_flag >= 5)
    config.plot.f2  = f2;
    config.plot.f3  = f3;
end

config.tracking = tracking;
config.trackingnr = 1; 
 
% Image parameters
Image.flag = image_flag; 
Image.sizes = [size_x, size_y, size_z]; 

% Information for initial curve generation
initial_information.center = c; 
initial_information.radius = r; 
initial_information.J      = J; 

% Surface parameters
Surface.flag = surface_flag; 
Surface.R = R;
Surface.r = r0; 
Surface.triangulation = triangulation;  
Surface.tol = surface_tol;
Surface.level_max = level_max;

end
