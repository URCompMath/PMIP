function Icb = my_rgb2cb(I)

I = (I)/255; 

chrom = [I(1),I(2),I(3)]; 
bright = norm(chrom); 

if(bright > 0)
    chrom = chrom /bright;
end

Icb = [chrom,bright]'; 
end