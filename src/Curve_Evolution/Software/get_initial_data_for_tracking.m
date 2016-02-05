function [Gamma, Omega] = get_initial_data_for_tracking(config, Image)

%% Image data
if(config.dimension == 2)
    switch Image.flag
        
        case 6 % Satellite Tracking
            filename ='../Input/2d_tracking/initial_data.mat';
            load(filename);
            
    end
% Here: Possibility to call other input data for other tracking examples

end
