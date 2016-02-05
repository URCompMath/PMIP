function Ihsv = my_rgb2hsv(I)

red = (I(1))/255;
green = (I(2))/255;
blue = (I(3))/255;

% red,green,blue in [0,1]
% compute h \in [0,2pi], s,v \in [0,1]

rgb_min = min([red,green,blue]);
rgb_max = max([red,green,blue]);

val = rgb_max;
delta = rgb_max - rgb_min;

if(rgb_max > 0)
    sat = delta / rgb_max;
    
    if(delta > 0)
        if(red == rgb_max)  % between yellow and magenta
            hue = (green - blue)/ delta;
        else
            if(green == rgb_max) % between cyan and yellow
                hue = 2+ (blue - red)/ delta;
            else % between magenta and cyan
                hue = 4+ (red - green)/delta;
            end
        end
    else
        % delta = 0 --> rgb_max = rgb_min --> red=green=blue --> gray color
        hue = 0;
        sat = 0;
    end
    
    % Transform to radians
    hue = hue * pi/3;
    
    % Shift negative angles by +2pi
    if(hue < 0)
        hue = hue + 2*pi;
    end
    
else
    % red=green=blue=0, sat=0, val is arbitrary, set to 0
    sat = 0;
    hue = 0;
end

Ihsv = [hue;sat;val]; 
end