% init  parameters

% PDE Flag
config.method.pde_flag = 1;  % 0: mean curvature flow (MCF) 1: image segmentation (default: image segmentation)

% Image Flag
config.image.flag = 1;   

% 1  Sphere MCF                      
% 2  Splitting Two Balls             
% 3  Merging to one ball            
% 4  Ball 2 Torus                   
% 5  Torus 2 Ball                 
% 6  Lung TCIA                    
% 7  UKR Lung Segm              
% 8  UKR Amdominal Part         
% 9  UKR Lung Splitting 

%% Default parameters
% Sizes
config.image.sizes = [5,5,5];
config.image.h = 5/50; 

% Plot rate
config.plot.nr = 1;

% Solve
config.solve.tol = 1e-1; 

% Smallest area before deletion
config.delete.area_tol = 4*pi*(0.1)^2;

% Refine, Coarsen Parameters
config.mesh.refinecoarsen.flag = 1; 
config.mesh.refine_factor = 10;    % refine simplex if area > factor * average area
config.mesh.factor_min = 1.5;
config.mesh.phi_max = 160*pi/180;  % refine simplex if one angle > phi_max
config.mesh.phi_min = 2*pi/180;    % delete simplex if one angle < phi_min
config.mesh.alpha_max = 20*pi/180;
config.mesh.area_desired = 1;      % for adaptive setting of refinement factor
config.mesh.improve_flag = 0; 

% Top change parameters
config.topchange.a = 0.05;
config.topchange.Ndetect = 10;
config.topchange.flag = 1; 
config.topchange.merge_angle = 150*pi/180; 
config.topchange.thr = [20,150,40]*pi/180;
config.topchange.countermin = 5;
config.topchange.a_factor = 10;

% Time step size control parameters
config.timestepsize.delta_t_init = 0.001;
config.timestepsize.flag = 1; 
config.timestepsize.factor_delta_t = 1;
config.timestepsize.decrease_factor = 10;
config.timestepsize.tol_DXn_min = 0.005;
config.timestepsize.tol_DXn_max = 0.1;


% Mesh regularization parameters
config.mesh.flag = 2;    % 1: not used
                         % 2: Induced Tang. Motion of Barrett, Garcke,
                         %    Nürnberg (strategy (i) (lambda_t = 0) or
                         %    (ii) (lambda_t = 1)
                         % 0: no mesh regularization
                         
config.mesh.sigma_t = 1;  % BGN
config.mesh.lambda_t = 1;
config.mesh.Niter = 8;    
config.mesh.change_tol = 0.1; 
config.mesh.delete_flag = 0;


% Nr of surfaces and max Nr of surfaces
Gamma.nr_surfaces = 1;
config.max_nr_surfaces = 100;

% Initial information, only needed for image segmentation (not for mean
% curvature flows)
initial_info.flag = 2;     % 1: sphere(s)/ellipsoid(s), 2: quader(s), cube(s), 3: 2 balls
initial_info.nr   = 1;     % number of initial surfaces
initial_info.data = [0,0,0,1,1,1]; % default data

% Region update parameters
config.width_band = 1;

% Multiphase flag
config.multiphase = 0; 


%% Image dependent parameters
switch config.image.flag
    case 1    % Mean Curvature Flow of a dumpbell with pinch off
        config.method.pde_flag = 0; % MCF
        config.initial_surface_filename = '../Input/initial_data_pinchoff';
                         
        config.image.sizes = [3,1.5,1.5];
        config.image.h = 6/120; 

        config.method.pde_flag = 0; 
        config.K = 300;

        config.method.sigma = 1;
        config.method.lambda = 0; 
        
        config.mesh.flag = 2;  
        config.mesh.refinecoarsen.flag = 0; 
        config.mesh.sigma_t  = 0.1; 
        config.mesh.lambda_t = 0;
        config.mesh.delete_flag = 0;
 
        config.topchange.flag = 1; 
        config.topchange.a = 0.025;  
        config.topchange.Ndetect = 10;
        config.topchange.thr = [30,150,40]*pi/180;
        config.topchange.countermin = 2;

        config.timestepsize.flag = 1; 
        config.timestepsize.tol_DXn_min = 0.003; 
        config.timestepsize.tol_DXn_max = 0.05;
        
        
     case 2   % Image Segmentation, Splitting Example, 1 Cylinder to 2 Balls
        config.initial_surface_filename = '../Input/initial_data_cylinder'; 
        
        config.K = 300;
        
        config.image.sizes = [2.1 2.1 2.1];
        config.image.h = 1/80; 
        
        initial_info.flag = 5;   % Initial Surface: Cylinder
        initial_info.data = [0,0,0,0.8,3.2,0]; % centered at origin, r=0.8, h=3.2
        
        config.method.sigma = 1;
        config.method.lambda = 100; 

        config.mesh.flag = 2; 
        config.mesh.refinecoarsen.flag = 0; 
        config.mesh.sigma_t  = 0.1; 
        config.mesh.lambda_t = 0;
        config.mesh.alpha_max = 30*pi/180;
        config.mesh.delete_flag = 0;
        
        config.topchange.flag = 1; 
        config.topchange.a = 0.025; 
        config.topchange.Ndetect = 10;
        config.topchange.thr = [30,150,40]*pi/180;
        config.topchange.countermin = 5;

        config.timestepsize.flag = 1; 
        config.timestepsize.tol_DXn_min = 0.003;
        config.timestepsize.tol_DXn_max = 0.04;
        
        
    case 3    % Image Segmentation, 2 balls merging
        config.initial_surface_filename = '../Input/initial_data_2balls_4merging'; 

        config.K = 200;
        
        config.image.sizes = [1.1 1.1 1.1];
        config.image.h = 1/80; 
        
        initial_info.flag = 3;    % Initial Surface: Two Balls
        initial_info.data = [0.5,0,0,0.4,0.4,0.4; -0.5,0,0,0.4,0.4,0.4;]; % Centered at +- [0.5,0,0], rx=ry=rz=0.4
        Gamma.nr_surfaces = 2;

        config.method.sigma = 2;
        config.method.lambda = 60;
        
        config.mesh.flag = 2;
        config.mesh.refinecoarsen.flag = 1; 
        config.mesh.refine_factor = 4;    
        config.mesh.sigma_t  = 10; 
        config.mesh.lambda_t = 0;
        config.mesh.alpha_max = 30*pi/180;
        config.mesh.area_desired = 0.001;
        config.mesh.factor_min = 1.5;
        config.mesh.phi_min = 2*pi/180;   
        config.mesh.phi_max = 170 * pi/180;
        config.mesh.delete_flag = 1;
        config.mesh.improve_flag = 1;
        
        config.topchange.flag = 1; 
        config.topchange.a = 0.0565/2;
        config.topchange.Ndetect = 8;
        
        config.timestepsize.flag = 1; 
        config.timestepsize.tol_DXn_min = 0.001; 
        config.timestepsize.tol_DXn_max = 0.02;
        
        
        
    case 4  % Image Segmentation, ball 2 torus (increase of genus example)
        config.initial_surface_filename = '../Input/initial_data_ellipsoid'; 
       
        config.method.pde_flag = 1; % Chan-Vese image segmentation
        config.K = 1000;
        
        config.image.sizes = [2.1 2.1 2.1];
        config.image.h = 1/80; 
        
        initial_info.flag = 1;   % Initial Surface: Sphere/Ellipsoid
        initial_info.data = [0,0,0,1.5,1.5,1.5]; % centered at origin, r=1.5
        
        config.method.sigma = 1;
        config.method.lambda = 60; 
        
        config.mesh.flag = 2; 
        config.mesh.refinecoarsen.flag = 1; 
        config.mesh.refine_factor = 4;    
        config.mesh.sigma_t  = 10; 
        config.mesh.lambda_t = 0;
        config.mesh.alpha_max = 50*pi/180;
        config.mesh.area_desired = 0.005;
        config.mesh.factor_min = 1.5;
        config.mesh.phi_min = 2*pi/180;   
        config.mesh.phi_max = 170 * pi/180;
        config.mesh.delete_flag = 1;
        
        config.topchange.flag = 1; 
        config.topchange.Ndetect = 8;        
        config.topchange.a = 0.0565;
        config.topchange.merge_angle = 178*pi/180;

        config.timestepsize.tol_DXn_min = 0.01;
        config.timestepsize.tol_DXn_max = 0.1;
        config.timestepsize.flag = 1; 

    
        
   case 5  % Torus to ball (decrease of genus example)
        config.initial_surface_filename = '../Input/initial_data_torus'; 

        config.K = 1000;
        
        config.image.sizes = [2.1 2.1 2.1];
        config.image.h = 1/80; 
        
        initial_info.flag = 4;    % Initial Surface: Torus
        initial_info.data = [0,0,0,0.9,0.5,0]; % centered at origin, R=0.9, r=0.5

        config.method.sigma = 1;
        config.method.lambda = 20; 
        
        config.mesh.flag = 2;  
        config.mesh.refinecoarsen.flag = 0; 
        config.mesh.refine_factor = 10; 
        config.mesh.sigma_t  = 0.1; 
        config.mesh.lambda_t = 0;
        config.mesh.alpha_max = 50*pi/180;
        config.mesh.area_desired = 0.005;
        config.mesh.factor_min = 1.5;
        config.mesh.phi_min = 2*pi/180;   
        config.mesh.phi_max = 170 * pi/180;
        config.mesh.delete_flag = 1;

        config.topchange.flag = 1; 
        config.topchange.Ndetect = 20;
        config.topchange.a = 0.025; 
        config.topchange.merge_angle = 120*pi/180; 
        
        config.timestepsize.flag = 1; 
        config.timestepsize.tol_DXn_min = 0.0005;
        config.timestepsize.tol_DXn_max = 0.01;
        
        
    case 6 % Lung, from TCIA, Cancer Imaging Archive
        config.initial_surface_filename = '../Input/initial_data_lung_2ellipsoid'; 
        
        config.method.pde_flag = 1; 
        config.K = 1500;
        
        config.image.sizes = [445 310 250];
        config.image.h = 1; 
        
        initial_info.flag = 1;   % Two balls centered at c1 and c2 with radius r1 and r2
        initial_info.data = [455/3,310/2,250/2,65,100,75;   455*2/3,310/2,250/2,65,100,75]; 
        Gamma.nr_surfaces = 2;

        config.cross.y_const = [80 150 220]'; % Draw cross sections
        config.cross.z_const = [80 150 200]'; 

        config.method.sigma = 10;
        config.method.lambda = 400; 
        
        config.solve.tol = 5000;

        config.mesh.flag = 2;
        config.mesh.refinecoarsen.flag = 1; 
        config.mesh.refine_factor = 4;    
        config.mesh.sigma_t  = 10; 
        config.mesh.lambda_t = 0;
        config.mesh.alpha_max = 50*pi/180;
        config.mesh.area_desired = 25; 
        config.mesh.factor_min = 1.5;
        config.mesh.phi_min = 2*pi/180;   
        config.mesh.phi_max = 170 * pi/180;
        config.mesh.delete_flag = 0;
        
        config.topchange.flag = 1; 
        config.topchange.a = 1;
        config.topchange.Ndetect = 10; 
        config.topchange.thr = [30,150,40]*pi/180;
        config.topchange.countermin = 5;
        
        config.timestepsize.flag = 1;
        config.timestepsize.delta_t_init = 5;
        config.timestepsize.tol_DXn_min = 0.5; 
        config.timestepsize.tol_DXn_max = 10; 
        

        
    case 7 % UKR, lung segmentation
        config.initial_surface_filename = '../Input/initial_data_ukr_lung'; 
        
        config.K = 1500;
        
        config.image.sizes = [128 128 141];
        config.image.h = 1; 
        
        initial_info.flag = 1;   % Two balls centered at c1 and c2 with radius r1 and r2
        initial_info.data = [32 64 90 10 10 30;  94 64 90 10 10 30]; 
        Gamma.nr_surfaces = 2;
        
        config.cross.y_const = [38 50 68]'; % Draw cross sections
        config.cross.z_const = [10 40 90]'; 

        config.method.sigma = 1;
        config.method.lambda = 20; 
        
        config.solve.tol = 5000;

        config.mesh.flag = 2;
        config.mesh.refinecoarsen.flag = 1; 
        config.mesh.refine_factor = 4;    
        config.mesh.sigma_t  = 10; 
        config.mesh.lambda_t = 0;
        config.mesh.alpha_max = 50*pi/180;
        config.mesh.area_desired = 3; 
        config.mesh.factor_min = 1.5;
        config.mesh.phi_min = 2*pi/180;   
        config.mesh.phi_max = 170 * pi/180;
        config.mesh.delete_flag = 0;
        
        config.topchange.flag = 0; 
        config.topchange.a = 1;
        config.topchange.Ndetect = 10; 
        config.topchange.thr = [30,150,40]*pi/180;
        config.topchange.countermin = 5;
        
        config.timestepsize.flag = 1;
        config.timestepsize.delta_t_init = 1;
        config.timestepsize.tol_DXn_min = 0.25*0.2; 
        config.timestepsize.tol_DXn_max = 0.25*8; 
        
        
  
    case 8 % UKR, Abdominal Region
        config.initial_surface_filename = '../Input/initial_data_ukr_adominal'; 
        
        config.multiphase = 1;
        config.K = 100;
        
        config.image.sizes = [128 128 80];
        config.image.h = 1; 
        
        initial_info.flag = 1;   % 3 ellipsoids centered at ci with radius ri, i=1,2,3 
        initial_info.data = [ 30 45 40 10 10 40; 100 45 40 10 10 40; 70 60 40 6 6 40];
        Gamma.nr_surfaces = 3;
                      
        config.cross.y_const = [52 64 76]'; % Draw cross sections
        config.cross.z_const = [30 50 70]';
        
        config.method.sigma = 1;
        config.method.lambda = 1000;
        
        config.solve.tol = 5000;
        
        config.mesh.flag = 2; 
        config.mesh.refinecoarsen.flag = 1;
        config.mesh.refine_factor = 4;
        config.mesh.sigma_t  = 10;
        config.mesh.lambda_t = 0;
        config.mesh.alpha_max = 120*pi/180; 
        config.mesh.area_desired = 3;
        config.mesh.factor_min = 1.5;
        config.mesh.phi_min = 2*pi/180;
        config.mesh.phi_max = 170 * pi/180;
        config.mesh.delete_flag = 1;
        
        config.topchange.flag = 0;
        config.topchange.a = 2;
        config.topchange.Ndetect = 8;
        config.topchange.thr = [30,150,40]*pi/180;
        config.topchange.countermin = 5;
        config.topchange.merge_angle = 170*pi/180;
        
        config.timestepsize.flag = 1;
        config.timestepsize.delta_t_init = 1;
        config.timestepsize.tol_DXn_min = 0.2;
        config.timestepsize.tol_DXn_max = 5;
        
        
    case 9 % UKR, Lung Segmentation, Splitting Example
        config.initial_surface_filename = '../Input/initial_data_ukr_lung_splitting'; 
        
        config.K = 1000;
        
        config.image.sizes = [128 128 280];
        config.image.h = 1; 
        
        initial_info.flag = 1;   % ellipsoid centered at ci with radius ri, i=1,2,3 
        initial_info.data = [ 64 64 120 30 10 10]; 
        Gamma.nr_surfaces = 1;
                           
        config.cross.y_const = [50 64 80]'; % Draw cross sections
        config.cross.z_const = [80 120 160]';
        
        config.method.sigma = 1;
        config.method.lambda = 20; 
        
        config.solve.tol = 2000;
        
        config.mesh.flag = 2; 
        config.mesh.refinecoarsen.flag = 1;
        config.mesh.refine_factor = 4;
        config.mesh.sigma_t  = 100;
        config.mesh.lambda_t = 0;
        config.mesh.alpha_max = 120*pi/180; 
        config.mesh.area_desired = 3;
        config.mesh.factor_min = 1.5;
        config.mesh.phi_min = 2*pi/180;
        config.mesh.phi_max = 170 * pi/180;
        config.mesh.delete_flag = 1;
        
        config.topchange.flag = 1;
        config.topchange.a = 2;
        config.topchange.Ndetect = 8;
        config.topchange.thr = [30,150,40]*pi/180;
        config.topchange.countermin = 5;
        config.topchange.merge_angle = 170*pi/180;
        
        config.timestepsize.flag = 1;
        config.timestepsize.delta_t_init = 0.2;
        config.timestepsize.tol_DXn_min = 0.1;
        config.timestepsize.tol_DXn_max = 2;

end

config.mesh.sigma_t_init = config.mesh.sigma_t;



