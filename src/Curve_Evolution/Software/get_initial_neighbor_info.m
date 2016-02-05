function Gamma_new = get_initial_neighbor_info(Gamma)

J = size(Gamma.X,1);
Gamma.neigh = zeros(J,2);

Nmain = size(Gamma.index_info,1);
Nsub  = size(Gamma.index_info,2);

for k=1:Nmain
    for l=1:Nsub
        i0 = Gamma.index_info(k,l);
        if(l<Nsub)
            iN = Gamma.index_info(k,l+1)-1;
        else
            if(k<Nmain)
                iN = Gamma.index_info(k+1,1)-1;
            else
                iN = J;
            end
        end
        
        % Assume closed curves, no boundary contact at the beginning
        % i=i0
        Gamma.neigh(i0,1)=iN;
        Gamma.neigh(i0,2)=i0+1;
        % i=iN
        Gamma.neigh(iN,1)=iN-1;
        Gamma.neigh(iN,2)=i0;
        % i=i0+1:iN-1
        for i=(i0+1):(iN-1)
            Gamma.neigh(i,1)=i-1;
            Gamma.neigh(i,2)=i+1;
        end
    end
end

Gamma_new = Gamma;
end