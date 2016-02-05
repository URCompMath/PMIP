function I_out = func_I(z,image_flag,Omega)

switch image_flag
   
 %     case 1 % not used here mcf
       
    case 2  % two balls
        c = 1.2;
        r = 0.8;
        if(norm(z-[-c;0;0])<=r || norm(z-[c;0;0])<=r)
            I_out = 0;
        else 
            I_out = 255; 
        end
     
        
    case 3  % merging
        c = 0; 
        r = 0.6; 
        if(norm(z-[c;0;0])<=r)
            I_out = 0;
        else 
            I_out = 255; 
        end       
       
         
    case 4  % torus
        r=0.4;
        R=1.2;

        if((sqrt(z(1)^2+z(2)^2)-R)^2+z(3)^2<=r^2)
            I_out = 0;
        else 
            I_out = 255;
        end
        
        
    case 5 % torus 2 ellipsoid
        r = [0.8,0.8,0.8];
        
        if((z(1)/r(1))^2 + (z(2)/r(2))^2 + (z(3)/r(3))^2 <= 1)
            I_out = 0;
        else 
            I_out = 255; 
        end 
        

    case 6 % Lung, Cancer Imaging Archive
        k = 310+1-round(z(2)); 
        l = round(z(1)); 
        m = round(z(3)); 
        
        if(k<= 0 || k>310 || l<=0 || l>445 || m<=0 || m>250)
            I_out = 160;
        else
            I_out = double(Omega.data(k,l,m));
        end
        
        
    case 7 % UKR, lung 
        k = 128+1-round(z(2)); 
        l = round(z(1)); 
        m = round(z(3)); 
        
        if(k<= 0 || k>128 || l<=0 || l>128 || m<=0 || m>141)
            I_out = 100; 
        else
            if(Omega.Omega_info(k,l,m) == 0)
                I_out = 100;
            else
                I_out = double(Omega.data(k,l,m)); 
            end
        end
 
        
    case 8 % UKR, abdominal part
        k = 128+1-round(z(2)); 
        l = round(z(1)); 
        m = round(z(3)); 
        
        if(k<= 0 || k>128 || l<=0 || l>128 || m<=0 || m>80)
            I_out = 70; 
        else
            if(Omega.Omega_info(k,l,m) == 0)
                I_out = 70;
            else
                I_out = double(Omega.data(k,l,m)); 
            end
        end
        
        
    case 9 % UKR, lung splitting
        k = 128+1-round(z(2)); 
        l = round(z(1)); 
        m = round(z(3)); 
        
        if(k<= 0 || k>128 || l<=0 || l>128 || m<=0 || m>280)
            I_out = 70; 
        else
            if(Omega.Omega_info(k,l,m) == 0)
                I_out = 70;
            else
                I_out = double(Omega.data(k,l,m)); 
            end
        end
end

end

