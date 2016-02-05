function mark = mark_simplices_for_deletion_collision(X,simplices,nodes)

N = size(simplices,1); 
mark = zeros(N,1); 

Nnodes = size(X,1); 
markX = zeros(Nnodes,1); 

for i=1:size(nodes,1)
    markX(nodes(i,1),1) = 1;
end


for i=1:N
    for j=1:3
        i0 = simplices{i,1}.nodes(j); 
        if(markX(i0,1) == 1)
            mark(i,1) = 1; 
        end
    end
end
Ndelete = sum(mark); 
fprintf('%d simplices marked for deletion due to collision\n', Ndelete); 

end
