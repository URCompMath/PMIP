function [config,Image,Surface,initial_information] = init_parameters_2d()
% PARAMETER FILE

% Image flag
image_flag = 1; 
% 1: ABC
% 2: Multiphase Gray (3 objects)
% 3: Triple Junctions and Color
% 4: Medical Image
% 5: Flowers
% 6: Satellite Tracking


%% Set the parameters

%% Default parameters
% (will be partly overwritten dependent on image_flag)

% Image Segmentation Method
im_seg_method = 1;   % 1: Chan Vese with parametric method incl. multiple
% phases, incl. triple junctions and topology changes
% (parameter im_seg_method can be used to implement other segmentation
% methods like Geodesic Active Contours here, currently only Chan-Vese)


% Lmax and Lmin, parameters for global refinement/ coarsening
Lmax =  5; 
Lmin =  2.5; 
Nmin = 10;

% Use different length tolerance for new curves (new small curves should
% not be immediately deleted, time for growing --> age_min time steps)
age_min = 5; 

% Width of narrow band where the regions are updated around Gamma
width_band = 1;

% Substeps if main step is divided
n_sub = 10;
k_top_last = -100;
kwait = 2;

% Grid size top change
amin = 1;
amax = 4;

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

% Band width for Prev in detect_top_change
a_Prev = 10;

% Setting of sigma
sigma_factor = 0.2;
sigma_compute = 10;
sigma_min = 30;
sigma_step = 5;


%% Image dependent parameters
switch(image_flag)
    case 1 % ABC
        % Image sizes
        size_x     = 225;
        size_y     = 125;
        % number of nodes (initial number)
        J= 15*32;
        % Number of timesteps
        K=600; 
        % Main time step size
        Delta_t = 0.1;
        % Parameter sigma (weight for curvature)
        sigma = 1; 
        sigma_compute = 0;
        sigma_min = 1;
        sigma_factor = 0.2;
        % Parameter lambda (weight for external force)
        lambda = 20; 
        % Parameter lambda partial
        lambda_partial = 0.8*sigma;
        % Parameter diffusion  - mu Laplace u + u = u_0
        param_diffusion = 5;
        % Center and radius of intial contour(s)
        c = zeros(2,2,1); 
        c(:,1,1) = [70; 125/2]; 
        c(:,2,1) = [150; 125/2]; 
        r = 25*ones(2,2,1);         
        % Check topological changes? 0: no, 1: yes
        top_check = 1;
        % Length tolerance
        length_tol = 2;
        % Plot number
        plot_nr = 10;
        plot_nr2 = 10;
        % Grid size top change
        amin = 1;
        amax = 4; 
        % Refine / Coarsening Parameters
        Lmax =  4;
        Lmin =  1; 
 
        
    case 2 % Multiphase Test Image: 3 objects
        % Image sizes
        size_x     = 300;
        size_y     = 300;
        % number of nodes (initial number)
        J=300;
        % Number of timesteps
        K=2600;
        % Main time step size
        Delta_t = 0.1;
        % Parameter sigma (weight for curvature term)
        sigma = 5; 
        sigma_min = 1;
        sigma_factor = 0.3;
        % Parameter lambda (weight for external force)
        lambda = 20; %5;
        % Parameter lambda partial
        lambda_partial = 0.8*sigma;
        % Parameter diffusion  - mu Laplace u + u = u_0
        param_diffusion = 5;
        % Center and radius of intial contour(s)
        c = zeros(2,1,3);
        c(:,1,1) = [size_x/3; size_y*2/3];
        c(:,1,2) = [2/3*size_x; size_y*2/3];
        c(:,1,3) = [size_x/2; size_y/3];
        
        r = 40*ones(2,1,3);
        
        % Check topological changes? 0: no, 1: yes
        top_check = 1;
        % Length tolerance
        length_tol = [2,5];
        % Plot number
        plot_nr = 10;
        plot_nr2 = 10;

        
    case 3 % Colored Test Image, 3 colored balls with triple junctions
        % Image sizes
        size_x     = 300;
        size_y     = 300;
        % number of nodes (initial number)
        J=90;
        % Number of timesteps
        K=300;
        % Main time step size
        Delta_t = 0.1;
        % Parameter sigma (weight for curvature term)
        sigma = 30;
        sigma_compute = 0; 
        % color method
        color_method = 1;  % 1: RGB (red-green-blue channels)
                           % 2: CB (Chromaticity Brightness)
                           % 3: HSV (Hue-Saturation-Value)
        % Parameter lambda (weight for external force)
        lambda = 5;
        % Parameter lambda partial 
        lambda_partial = 0.8*sigma;
        % Parameter diffusion  - mu Laplace u + u = u_0
        param_diffusion = 5;
        % Center and radius of intial contour(s)
        c = zeros(2,1,3);
        c(:,1,1) = [100;160];     % Gamma 1
        c(:,1,2) = [150;130];     % Gamma 2
        c(:,1,3) = [200,160];     % Gamma 3  
        
        r = 25*ones(2,1,3);
        
        % Check topological changes? 0: no, 1: yes
        top_check = 1;
        % Length tolerance
        length_tol = 2;
        % Plot number
        plot_nr = 10; 
        plot_nr2 = 10;
        
        
    case 4  % Medical image
        % Image sizes
        size_x     = 275;
        size_y     = 275;
        % number of nodes (initial number)
        J=25*32; 
        % Number of timesteps
        K=700;
        % Main time step size
        Delta_t = 0.02;
        % Parameter sigma (weight of curvature term)
        sigma = 20; 
        % Parameters for automatic setting of sigma
        sigma_compute = 50;
        sigma_min  = 5; 
        sigma_factor = 0.05; 
        sigma_step = 5;
        % Parameter lambda (weight for external force)
        lambda = 400;
        % Parameter lambda partial 
        lambda_partial = 0.8*sigma;
        % Parameter diffusion  - mu Laplace u + u = u_0  
        param_diffusion = 1;
        % Center and radius of intial contour(s)
        c = zeros(2,1,3);   
        c(:,1,4) = [17;130];
        c(:,1,1) = [55;130];
        c(:,1,2) = [100; 130];
        c(:,1,3) = [150; 130];
        c(:,2,4) = [200; 130];
        c(:,2,1) = [240; 130];
        
        c(:,2,2) = [30;40];
        c(:,2,3) = [80;40];
        c(:,3,4) = [135; 40];
        c(:,3,1) = [190; 40];
        c(:,3,3) = [250; 40];
        
        c(:,4,1) = [30;230];
        c(:,3,2) = [80;230];
        c(:,4,4) = [135;230];
        c(:,4,2) = [190; 230];
        c(:,4,3) = [250; 230];
        
        r = 15*ones(2,4,4);   
        
        % Check topological changes? 0: no, 1: yes
        top_check = 1;
        % Length tolerance
        length_tol = [1,5];
        % Plot number
        plot_nr = 10;
        plot_nr2 = 50;
        % Grid size top change
        amin = 1;
        amax = 2;
        % Refine / coarsening parameters
        Lmin = 1.5; 
        Lmax = 4;
        % Set age min 
        age_min = 15; 
    
    
    case 5 % Flowers Berkley
        % Image sizes
        size_x     = 481;
        size_y     = 321;
        % number of nodes (initial number)
        J=15*32;
        % Number of timesteps
        K=650;
        % Main time step size
        Delta_t = 0.01;
        % Parameter sigma (weight for curvature)
        sigma = 125;
        sigma_min = 20;
        sigma_factor = 0.15; 
        sigma_compute = 10;        
        % color method
        color_method = 2;  % 1: RGB (red-green-blue channels)
                           % 2: CB (Chromaticity Brightness)
                           % 3: HSV (Hue-Saturation-Value)
        
        % Parameter lambda (weight for external force)
        % CB method lambda = [lambda_c, lambda_b]
        lambda = [180, 40];
        % Parameter lambda partial 
        lambda_partial = 0.4*sigma;
        % Parameter diffusion  - mu Laplace u + u = u_0
        param_diffusion = 1;
        % Center and radius of intial contour(s)    
        ay = size_y/6;
        ax = size_x/10;

        r = 45*ones(2,3,5);
        
        c = zeros(2,3,5);
        c(:,1,2) = [ax;1*ay];
        c(:,1,4) = [3*ax;1*ay];
        c(:,1,1) = [5*ax;1*ay];
        c(:,1,5) = [7*ax;1*ay];
        c(:,3,3) = [9*ax;1*ay];
        
        c(:,2,5) = [ax;3*ay];
        c(:,1,3) = [3*ax;3*ay];
        c(:,3,2) = [5*ax;3*ay];
        c(:,2,3) = [7*ax;3*ay];
        c(:,2,4) = [9*ax;3*ay];
        
        c(:,3,4) = [ax;5*ay];
        c(:,2,1) = [3*ax;5*ay];
        c(:,3,5) = [5*ax;5*ay];
        c(:,3,1) = [7*ax;5*ay];
        c(:,2,2) = [9*ax;5*ay];
        
        for i=1:3
            for j=1:5
                c(2,i,j)=size_y - c(2,i,j); 
            end
        end

        % Check topological changes? 0: no, 1: yes
        top_check = 1;
        % Length tolerance
        length_tol = [2,9];
        % Plot number
        plot_nr = 10;
        plot_nr2 = 10; 
        % Refine / Coarsening Parameters
        Lmax =  6; 
        Lmin =  2.5; 
        % Grid size top change
        amin = 1;
        amax = 3;
        

             
    case 6   % CamSim Tracking
        % Image size
        size_x = 400;
        size_y = 400;
        % number of nodes (initial number)
        J=200; % overwritten when initial curve is loaded during tracking!
        % Number of timesteps per image
        K=20;
        % Main time step size
        Delta_t = 0.02;
        % Parameter sigma (weight for curvature)
        sigma = 50;
        sigma_compute = 0;
        % Parameter lambda (weight for external force)
        lambda = 80;
        % Parameter lambda partial 
        lambda_partial = 0.8*sigma;
        % Parameter diffusion  - mu Laplace u + u = u_0
        param_diffusion = 1;
        % Center and radius of intial contour(s)
        c = zeros(2,2,2);  % is overwritten when initial curve is loaded during tracking!
        c(:,1,1)=[230;100];
        c(:,2,1)=[230;240];
        c(:,1,2)=[150;200];
        c(:,2,2)=[280;200];
 
        r = 25*ones(2,2,2);
        % Check topological changes? 0: no, 1: yes
        top_check = 1;
        % Length tolerance
        length_tol = [2,9];
        % Plot number
        plot_nr = 5;
        plot_nr2 = 10;
        % Grid size top change
        amin = 1;
        amax = 4; 
        Nmin = 20;
        % Refine / Coarsening Parameters
        Lmax =  6;
        Lmin =  2;      
        % Tracking number
        tracking = 6;
end


f1=figure;
set(f1, 'Position', [136   373   560   420]);



%% Write parameters to structure

% Configuration parameters
config.dimension = 2;
config.K = K;
config.timestep = Delta_t;

config.method.flag = im_seg_method;
config.method.sigma = sigma;
config.method.sigma_factor = sigma_factor;
config.method.sigma_compute = sigma_compute;
config.method.sigma_min = sigma_min;
config.method.sigma_step = sigma_step;
config.method.lambda = lambda;
config.method.lambda_partial = lambda_partial;
config.method.mu = mu; 

config.method.color = color_method;
config.param_diffusion = param_diffusion;

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

config.tracking = tracking;
config.trackingnr = 1;


% Image parameters
Image.flag = image_flag;
Image.sizes = [size_x, size_y];

% Information for initial curve generation
initial_information.center = c;
initial_information.radius = r;
initial_information.J      = J;

% Surface parameters
Surface = struct([]);  % empty surface, as dimension = 2;

end

