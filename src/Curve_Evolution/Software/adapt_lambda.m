function config_new = adapt_lambda(config,Omega,color_flag)

lambda = config.method.lambda;



if(size(lambda,1) ~= size(Omega.coeffs,1))
    lambda = ones(size(Omega.coeffs,1),1)* lambda(1,:);
end

if(color_flag)
    switch config.method.color
        case 1 % RGB
            % lambda = [lambda_red,lambda_green,lambda_blue]
            if(size(lambda,2) ~= 3) % Adaptation if one lambda for all channels
                lambda = lambda(:,1)*ones(1,3);
            end
        case 2 % CB
            % lambda = [lambda_c, lambda_b]
            if(size(lambda,2) ~= 2) % Adaptation if one lambda for chrom & brightness
                lambda = lambda(:,1)*ones(1,2);
            end
        case 3 % HSV
            % lambda = [lambda_h, lambda_s, lambda_v]
            if(size(lambda,2) ~= 3) % Adaptation if one lambda for h,s and v component
                lambda = lambda(:,1)*ones(1,3);
            end
    end
end

config_new = config;
config_new.method.lambda = lambda;
end
    
    