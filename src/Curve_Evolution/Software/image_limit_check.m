function Gamma_new = image_limit_check(Gamma,Image)
size_x = Image.sizes(1);
size_y = Image.sizes(2);

for i=1:size(Gamma.X,1)
    if(Gamma.X(i,1)<0.5)
        Gamma.X(i,1)=0.5;
    else
        if(Gamma.X(i,1)>size_x+0.5)
            Gamma.X(i,1)=size_x+0.5;
        end
    end
    if(Gamma.X(i,2)<0.5)
        Gamma.X(i,2)=0.5;
    else
        if(Gamma.X(i,2)>size_y+0.5)
            Gamma.X(i,2)=size_y+0.5;
        end
    end
end

Gamma_new = Gamma;
end