function Irgb = my_cb2rgb(I)

Irgb = [I(1),I(2),I(3)]'; 
Irgb = Irgb * I(4) * 255; 
end