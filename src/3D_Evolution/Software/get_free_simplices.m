function [simp_free, v_simp_free, index_simp_free] = get_free_simplices(Gamma)

simp_free = zeros(size(Gamma.simplices,1),3);   % 1st col: index of simplex, 2nd col: local index of free edge, 3rd column: group index
Nsimp_free = 0;
for i=1:size(Gamma.simplices,1)
    neigh = Gamma.simplices{i,1}.neigh;
    Nfree = 0;
    local_index = -1;
    for j=1:3
        if(neigh(j)==-1)
            Nfree = Nfree + 1;
            local_index = j;
        end
    end
    
    if(Nfree == 1 )
        Nsimp_free = Nsimp_free + 1;
        simp_free(Nsimp_free,:) = [i, local_index,-1];
    end
end
simp_free = simp_free(1:Nsimp_free,:);
v_simp_free = zeros(Nsimp_free,3);
index_simp_free = zeros(Nsimp_free,1);

fprintf('Free simplices:\n');
for i=1:Nsimp_free
    nodes = Gamma.simplices{simp_free(i,1),1}.nodes;
    X1 = Gamma.X(nodes(1),:);
    X2 = Gamma.X(nodes(2),:);
    X3 = Gamma.X(nodes(3),:);
    
    cross_j = cross(X1-X3,X2-X3);
    v  = cross_j/norm(cross_j);
    v_simp_free(i,:) = v;
    index_simp_free(i,1) = Gamma.simplices{simp_free(i,1),1}.index(2);
    
    fprintf('simp %d, local index  %d, neigh(%d)=%d, surface index  %d, nu = %2.3f %2.3f %2.3f\n', simp_free(i,1), simp_free(i,2), simp_free(i,2), Gamma.simplices{simp_free(i,1),1}.neigh(simp_free(i,2)), Gamma.simplices{simp_free(i,1),1}.index(2), v(1),v(2),v(3));
end
end