

function GFFDtest(h)


if h ==1
    disp("Attractive-Repulsive kernel")
    X = 2;
    T = 0.8;
    m = 8;
    k = 10;
    xs = linspace(-X,X,((2*m)+1));
    Pinit = arrayfun(@Pex31,xs)'; 
    GFFDSolve(@H0,@Wfunct,@V0,X,m,T,k,Pinit, @true31);

elseif h ==2
    disp("Nonlinear Diffusion with nonlocal attraction kernel")
    X = 3;
    T = 0.5;
    m = 15;
    k = 40;
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
    y = abX^2 *(1/(2 - log(abX)));
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
v = 1.48;
m = 3;
y = (v/m)*(x.^m);
end

function y = true32(x)% ---------- not given in paper
y = 0;
end

%----------------------------------------------------------GFFDSolve--------------------------------
function GFFDSolve(H,Wfun,V,X,m,T,k,P0, Soln)
time = linspace(0,T,k);
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
figure(1)
spaced =linspace(-X,X, M+49);
trueSoln = arrayfun(Soln,spaced);
tiledlayout(1,3)

nexttile
plot( space, P(:,1),'-o');
title("Initial Condition");


nexttile
plot(space, P(:,l-1), '-o');
title("Approximation");

nexttile
plot( spaced, trueSoln);
title("True Solution");


end

%---------------------------------------------------------FVF---------------------------------------------------

function Pnew = FVF(H,Wfun,V,X,m,T,k,pr)
M = 2*m;
space = linspace(-X,X,(M+1));
dt = T/k;
dx = 2*X/(M);
coeff = dt/dx;
%P = [P0, zeros(M+1,k)];
tester = 0;
theta = 2;
%P(:, 1) = P0;
%Evaluating V ahead to save computation
Vs = V(space)';
Z = Vs + H(pr)';
%create U
    %Need z_j = dx \sum_i W_{j-i}p_i  +H'(p_j) + V(x_j)   
    for j = 1:M+1
        con = space(j) - space;
        Wcon = arrayfun(Wfun,con);
        s = Wcon.*(pr');
        %disp(size(Wcon));
        wc = sum(s); 
        Z(j) = Z(j) + dx*wc;
    end
    Upos = zeros(M+1,1);
    Uneg = zeros(M+1,1);
    for j = 1:M+1
        if j ==M+1
            U = 0-Z(j);
        else
            U = Z(j+1) - Z(j);
        end
        if U>0
            Upos(j) = (1/dx)*U;
            Uneg(j) = 0;
        else 
            Upos(j) = 0;
            Uneg(j) = (1/dx)*U;
        end
    end

    MatA = (1/2)*eye(M+1) + coeff*(diag(Upos(1:M),-1)-diag(Upos));
    MatB = (1/2)*eye(M+1) - coeff*(diag(Uneg) + diag(Uneg(2:end),1));

    %Create (px) 
    px = (1/(2*dx))*(diag(ones(M,1),1)-diag(ones(M,1),-1))*pr;
    Pw = (pr - (dx/2)*px); %p^W
    Pe = (pr + (dx/2)*px); %p^E

    %if anything becomes negative change that entry according to minmod:
    for cn = 1:M+1
        if Pw(cn) <0 || Pe(cn)< 0
            tester = tester -1;
            if cn==1
                x = (theta/dx)*(pr(cn+1)-pr(cn));
                y = (pr(cn+1))/(2*dx);
                z = (theta/dx)*(pr(cn));
            elseif cn==M+1
                x = (theta/dx)*(-pr(cn));
                y = (0-pr(cn-1))/(2*dx);
                z = (theta/dx)*(pr(cn)-pr(cn-1));
            else
                x = (theta/dx)*(pr(cn+1)-pr(cn));
                y = (pr(cn+1)-pr(cn-1))/(2*dx);
                z = (theta/dx)*(pr(cn)-pr(cn-1));
            end
            Pw(cn) = (pr(cn) - (dx/2)*minmod(x,y,z));
            Pe(cn) = (pr(cn) + (dx/2)*minmod(x,y,z));
        end
       % if Pe(cn) <0
       %     tester = tester -1;
       %     if cn==1
       %         x = (theta/dx)*(pr(cn+1)-pr(cn));
       %         y = (pr(cn+1))/(2*dx);
       %         z = (theta/dx)*(pr(cn));
       %     elseif cn==M+1
       %         x = (theta/dx)*(-pr(cn));
       %         y = (0-pr(cn-1))/(2*dx);
       %         z = (theta/dx)*(pr(cn)-pr(cn-1));
       %     else
       %         x = (theta/dx)*(pr(cn+1)-pr(cn));
       %         y = (pr(cn+1)-pr(cn-1))/(2*dx);
       %         z = (theta/dx)*(pr(cn)-pr(cn-1));
       %     end
       %     Pe(cn) = (pr(cn) + (dx/2)*minmod(x,y,z));
       % end
    end

    pnew = MatA *Pe + MatB*Pw; 
    if nonneg(pnew)
        Pnew = pnew; %[0; pnew(2:M);0];  
    else 
        disp("error: Nonnegativity encountered");
        Pnew = 0;
    end
end

function w = minmod(x,y,z)
hold = sign([x,y,z]);
if hold == [1,1,1]
    w = min([x,y,z]);
elseif hold == [-1,-1,-1]
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