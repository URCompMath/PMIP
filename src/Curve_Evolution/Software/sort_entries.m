function v_out=sort_entries(v,value)
% assume v: 1xN vector

N=size(v,2);


mark = zeros(1,N);
for i=1:N
    if(abs(v(1,i)-value)<1e-12)
        mark(i)=1;
    end
end

w=v;
w(mark==1)=[];



v_out = [w, value*ones(1,size(v,2)-size(w,2))];
end