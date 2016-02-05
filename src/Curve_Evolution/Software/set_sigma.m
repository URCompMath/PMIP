function config_new = set_sigma(config,Eint,Eext,factor)

% internal energy Eint = sigma * lengths
% external energy Eext

% objective: choose sigma such that Eint = factor*Eext

Eint = Eint / config.method.sigma; 

if(Eint > 0)
    sigma_new = factor * Eext / Eint;
    sigma_new = max([sigma_new, config.method.sigma_min ]);
    fprintf('Change sigma from %3.2f to ', config.method.sigma); 
    config.method.sigma = sigma_new; 
    fprintf('%3.2f\n\n', config.method.sigma);
end

config_new = config; 
end


