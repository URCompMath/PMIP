function I_out = func_I(z,Image,Surface)

dim = max(size(z));
I = Image.data;

switch dim
    case 2
        Nx = size(I,2);
        Ny = size(I,1);
        
        k = min([Ny, max([round(Ny+1-z(2)),1])]);
        l = min([Nx, max([round(z(1))    , 1])]);
        
        if(size(size(I),2)==3)  % vector valued image (e.g. RGB colored)
            I_out = reshape(double(I(k,l,:)),size(I,3),1); 
            
        else
            I_out = double(I(k,l));     % scalar image (e.g. gray value)
        end
        
    case 3
          [cp,i_simp] = get_closest_point_and_simplex(z,Surface);
          I_out = Image.data(i_simp,:)'; 
end

end
