% init file
                      
switch dimension
    case 2
        [config,Image,Surface,initial_information] = init_parameters_2d();

    case 3
        [config,Image,Surface,initial_information] = init_parameters_geodesic();
    
    otherwise
        fprintf('Invalid dimension flag. Set dimension flag to 2.\n'); 
        dimension = 2; 
        [config,Image,Surface,initial_information] = init_parameters_2d();
end
