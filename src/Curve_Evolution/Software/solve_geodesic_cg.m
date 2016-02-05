function X=solve_geodesic_cg(S,b,config,repeat)

% Set tolerance
if(repeat == 0)
    tol = config.tol.std;
else
    tol = config.tol.extra;
end


kmax = 200;

% Precondion Matrix
C = diag(diag(S)); 

% initialize vectors
X = zeros(size(b));
r = - b + S*X; 
h = C*r;
d = -h; 

k=1;
res = r'*h; 

% Method of conjugate gradients
while(k<kmax && res > tol)
    Ad = S*d;
    alpha   = r'*h/(d'*Ad);
    X       = X + alpha*d;
    r       = r + alpha*Ad;
    h       = C*r;

    res_old = res;
    res     = r'*h;
    beta    = res/res_old;
    d       = -h + beta*d;
    
    k = k+1;
    
end
if(isnan(res))
    fprintf('Caution: res is NaN. Pausing computation.\n');
    pause
end
end