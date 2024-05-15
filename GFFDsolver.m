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
    Zpos = Vs + H(p);
    Zneg = Vs + H(p);
%Create px 
    for j = 2:M
        x = 2*(p(j+1)-p(j));
        y = (p(j+1)-p(j-1))/2;
        z = 2*(p(j)-p(j-1));
        px(j) = minmod(x,y,z);
    end
%create U
%Need z_j = dx \sum_i W_{j-i}p_i  +H'(p_j) + V(x_j)
    
    for j = 2:M
        wc = sum(wconv(j-1:(j+m-1)));
        Z = Zpos(j-1) + dx* wc;
        if Z>0
            Zpos(j-1) = Z;
            Zneg(j-1)=0;
        else 
            Zpos(j-1) = 0;
            Zneg(j-1)= Z;
        end
    end
    Up = diag(Zpos(2:M+1)-Zneg(1:M)) - diag(Zpos(2:M+1),-1) + diag(Zneg(1:M),1); 
    Upx = (-1)*diag(Zpos(2:M+1)-Zneg(1:M)) + diag(Zpos(2:M+1),-1) + diag(Zneg(1:M),1);
%Solve for next time step
    B = (dt/2) *Upx;
    A = eye(M) - coeff*Up;
    pnew = A*p - B*px;
    P(:,i+1) = pnew;  
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
