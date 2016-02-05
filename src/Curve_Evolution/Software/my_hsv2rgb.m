function Irgb = my_hsv2rgb(I)

H = I(1);
S = I(2);
V = I(3);


if(S==0)
    Irgb = [V;V;V];
else
    h = floor(H/(pi/3));
    f = H/(pi/3)-h;
    
    p = V*(1-S);
    q = V*(1-S*f);
    t = V*(1-S*(1-f));
    
    switch h
        case 0
            Irgb = [V;t;p];
        case 1
            Irgb = [q;V;p];
        case 2
            Irgb = [p;V;t];
        case 3
            Irgb = [p;q;V];
        case 4
            Irgb = [t;p;V];
        case 5
            Irgb = [V;p;q];
        case 6
            Irgb = [V;t;p];
    end
end

Irgb = Irgb * 255;
end