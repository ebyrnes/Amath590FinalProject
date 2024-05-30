

function GFFDtest_1(h)

if h ==1
    disp("Attractive-Repulsive kernel")
    X = 2;
    T = 3500;
    m = 10;
    k = 2800000;
    xs = linspace(-X,X,((2*m)+1));
    dx = 4/((2*m)+1);
    Pinit = zeros(size(xs)); %arrayfun(@Pex31,xs); 
    syms z
    Finit = (1/sqrt(2*pi)) *(exp(-(z^2)/2));
    for j = 1:((2*m)+1)
        xN = xs(j) - (dx/2);
       xP = xs(j) + (dx/2);
        Pinit(j) = (1.04767/dx)*int(Finit, xN, xP);% regularizing to make the mass 1 under the curve
    end
    
    GFFDSolve(@H0,@Wfunct,@V0,X,m,T,k,Pinit', @true31);

elseif h ==2
    disp("Nonlinear Diffusion with nonlocal attraction kernel")
    X = 3;
    T = 300;
    m = 15;
    k = 80000;
    xs = linspace(-X,X,((2*m)+1));
    Pinit = arrayfun(@Pex32,xs)'; 
    GFFDSolve(@H32,@W32,@V0,X,m,T,k,Pinit, @true32);
end

end

function y = Wfunct(x)
if x ==0
    y = 0;
else
    abX = abs(x);
    y = x^2 *(1/2) - log(abX);
end
end

function y = Pex31(x)
y = (1/sqrt(2*pi)) *(exp(-(x^2)/2));
end

function y = true31(x)
if abs(x) < sqrt(2)
    y = (1/pi) * (sqrt(2- (x^2)));
else
    y = 0;
end
end

function y = H0(j)
y = zeros(size(j));
end

function y = V0(j)
y = zeros(size(j));
end

function y = W0(j)
y = zeros(size(j));
end

function y = W32(x)
ep = (x^2)/(2); %2*theta in denom but theta = 1
y = -(1/sqrt(2*pi))*(exp(-ep));
end

function y = Pex32(x)
pos = exp(-(1/2)* (x+3)^2);
neg = exp(-(1/2)* (x-3)^2);
y = (1/(sqrt(8*pi))) *(pos +neg);
end

function y = H32(x)
%v = 1.48 & m = 3
y = (1.48)*(x.^(2));
end

%function y = true32(x)% ---------- not given in paper
%y = 0;
%end

%----------------------------------------------------------GFFDSolve--------------------------------
function GFFDSolve(H,Wfun,V,X,m,T,k,P0, Soln)
%time = linspace(0,T,k);
dt = T/k;
l =k+2;
M = 2*m;
space = linspace(-X,X,M+1);
P = [P0, zeros(M+1,k)];

for i = 1:k
    U = P(:,i);
    %first stage
    y1 = FVF(H,Wfun,V,X,m,T,k,U);
    if y1 == 0
        l =i;
        break
    end
    U1 = U + (dt)*y1;
    %second stage
    y2 =FVF(H,Wfun,V,X,m,T,k,U1);
    if y2==0
        l = i;
        break
    end
    U2 = (3/4)*U + (1/4)*U1 + (dt/4)*y2;
    %third stage
    y3 = FVF(H,Wfun,V,X,m,T,k,U2);
    if y3==0
        l = i;
        break

    end
    Un = (1/3)*U + (2/3)*U2 + (2/3)*dt* y3;
    P(:,i+1) = Un;
end

%disp("Solution over time/epochs");
%disp(P(:,1:k+1));
disp("Plotting")
disp(l)
%Plotting!

%spaced =linspace(-X,X, M+49);
disp("Final Solution is")
trueSoln = arrayfun(@true31,space);
disp(trueSoln);

disp("The error is: ")
err = P(:,l-1) - trueSoln';
disp(norm(err,2))

disp("The relative error is: ")
disp((norm(err,2)/norm(trueSoln,2)));

%figure(1)
%plot( space, P(:,1),'-o');
%title("Initial Condition");

%plot(space, err, '-o');
%title("Error Plot")

figure(1)
tiledlayout(1,2)

nexttile
plot(space, P(:,l-1), '-o');
title("Approximation");

nexttile
plot( space, trueSoln);
title("True Solution");
end

%---------------------------------------------------------FVF---------------------------------------------------

function Pnew = FVF(H,Wfun,V,X,m,T,k,pr)
M = 2*m;
space = linspace(-X,X,(M+1));
dt = T/k;
dx = 2*X/(M);
X1 = X+dx;
coeff = dt/dx;
%P = [P0, zeros(M+1,k)];
theta = 2;
%P(:, 1) = P0;
%Evaluating things ahead to save computation
prcn0 = pr(1); %ghostCell(H,V,space(1),dx,pr(1)); 
prcn1 = pr(end); %ghostCell(H,V,space(M+1),-dx,pr(M+1)); 
%Vs = V([space, X1])';
%Z = Vs + H([pr; prcn1])';
Z = zeros(M+2,1);
%create U
    %Need z_j = dx \sum_i W_{j-i}p_i  +H'(p_j) + V(x_j)   
    for j = 1:M+2
        if j ==M+2
            con = X1 - [space, X1];
        else
            con = space(j) - [space, X1];
        end
        Wcon = zeros(size(con)); 
        for q = 1:length(con) %Approximation of the integral
            W = (dx/10)*(Wfun(con(q)+(dx/2))+ 2*Wfun(con(q)+(dx/4))+ 2*Wfun(con(q)) + 2*Wfun(con(q)-(dx/4)) + Wfun(con(q)-(dx/2)));
            if q == M+2
                Wcon(q) = W * prcn1;
            else
                Wcon(q) = W * pr(q);
            end
        end
        %s = Wcon.*([pr; prcn1]');
        %disp(size(Wcon)); 
        wc = sum(Wcon); 
        Z(j) = Z(j) + dx*wc;
    end
    Upos = zeros(M+1,1);
    Uneg = zeros(M+1,1);

    for j = 1:M+1
            U = Z(j+1) - Z(j);
        if dt > dx/(2*abs(U))
            disp("CFL condition broken!")
            disp(U)
            break
        end
        if U>0
            Uneg(j) = -(1/dx)*U;
            Upos(j) = 0;
        else 
            Uneg(j) = 0;
            Upos(j) = -(1/dx)*U;
        end
    end

    MatAN = coeff*(diag(Upos(1:M),-1)); % for comp of F_j-1/2
    MatBN =  coeff*(diag(Uneg)); %changed 1:M -> 2:end and vice versa below
    
    MatAP = coeff*(diag(Upos)); % for comp of F_j+1/2
    MatBP = coeff*(diag(Uneg(2:end),1));

    %Create (px) 
   
    px  = (1/(2*dx))*(diag(ones(M+2,1),1)-diag(ones(M+2,1),-1))*[prcn0;pr;prcn1];
    pxC = px(2:M+2);

    Pw = pr - ((dx/2)*pxC); %p^W
    Pe = pr + ((dx/2)*pxC); %p^E
   
    %if anything becomes negative change that entry according to minmod:
    for cn = 1:M+1
        if Pw(cn) <0 || Pe(cn)< 0
            if cn==1
                x = (theta/dx)*(pr(cn+1)-pr(cn));
                y = (pr(cn+1)-prcn0)/(2*dx);
                z = (theta/dx)*(pr(cn)-prcn0);
            elseif cn==M+1
                x = (theta/dx)*(prcn1-pr(cn));
                y = (prcn1-pr(cn-1))/(2*dx);
                z = (theta/dx)*(pr(cn)-pr(cn-1));
            else
                x = (theta/dx)*(pr(cn+1)-pr(cn));
                y = (pr(cn+1)-pr(cn-1))/(2*dx);
                z = (theta/dx)*(pr(cn)-pr(cn-1));
            end
            Pw(cn) = (pr(cn) - (dx/2)*mnmd(x,y,z)); 
            Pe(cn) = (pr(cn) + (dx/2)*mnmd(x,y,z));
        end
    end
    FPl = MatAP * Pe + MatBP*Pw;
    FNg = MatAN * Pe + MatBN*Pw;
    FPl(end) = 0;
    FNg(1) = 0;
    Pnew = -1* (FPl - FNg); 
    %if nonneg(pnew)
    %    Pnew = pnew; %[0; pnew(2:M);0];  
    %else 
    %    disp("error: Nonnegativity encountered");
    %    Pnew = 0;
    %end
end
%Code from Ellie Byrnes


function [rhoOut] = ghostCell(Hp,V,x,deltax,rho)
target = Hp(rho) + V(x) - V(x-deltax);
targetfunc = @(r) Hp(r) - target;
rhoOut = fzero(targetfunc,rho);
end


function w = mnmd(x,y,z)
hold = sum(sign([x,y,z]));
if hold == 3
    w = min([x,y,z]);
elseif hold == -3
    w = max([x,y,z]);
else
    w = 0;
end
end

function t = nonneg(r)
t = true;
for i = 1:length(r)
    if r(i) <0
        t = false;
    end
end
end
