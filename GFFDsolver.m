%One-Dimensional Scheme

function GFFDsolver(H,W,V,X,m,T,k,P0)
time = linspace(0,T,k);
M = 2*m;
space = linspace(-X,X,M+1);
dt = T/k;
dx = 2*X/(M);
coeff = dt/dx;
P = zeros(M+1,(k+1));
P(:, 1) = P0;
%Evaluating functions ---------------------------- In vectorized notation rn!
Vs = V(space);
conv = linspace((-M*dx),(M*dx), M+1);
Wconv = W(conv);

%Forward Euler scheme to solve in time
for i = 1:k
    p = P(:, i);
    px = zeros(M-1,1);
    U = zeros(M-1,M-1);
    Z = zeros(M-1,1);
%Create px 
    for j = 2:M
        x = 2*(p(j+1)-p(j));
        y = (p(j+1)-p(j-1))/2;
        z = 2*(p(j)-p(j-1));
        px(j-1) = minmod(x,y,z);
    end
%create U
%Need z_j = dx \sum_i W_{j-i}p_i  +H'(p_j) + V(x_j)
    
    for j = 2:M
        wc = sum(wconv(j-1:(j+m-1)));
        Z(j-1) = wc;
    end
    Z = dx*Z + Vs + H(p);
    U = diag(Z) + diag(Z(),-1) + diag(Z(),1);%need to check shifts!-----------
%Solve for next time step
    B = coeff*U;
    A = eye(M) - B;
    pnew = A*p - B*px;
    P(:,+1) = pnew;  
end



end

%helper function for evaluating px above
function w = minmod(x,y,z)
hold = sign([x,y,z]);
if hold == [1,1,1]
    w = min(x,y,z);
elseif hold == [-1,-1,-1]
    w = max(x,y,z);
else
    w = 0;
end


end