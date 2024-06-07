function [xs,Ts,entropy2,entropy3,rho2,rho3] = test1DEntropy(Tmax)
h = 1;
if h == 1
    X = 4;
    m = 3;
    nu = 1.48;
    sigma = 1;
    theta=2;
    

    deltax = 0.1;
    deltat = deltax/10;
    k = ceil(Tmax/deltat);
    xs = -X:deltax:X;
    H = @(r) Hfunct2(r,nu,m);
    Wconv = Wfunct2(xs'-xs);
    rhoinit2 = zeros(size(xs));
    rhoinit3 = zeros(size(xs));
    rho0two = @(r) rho02(r,2);
    rho0three = @(r) rho02(r,3);
    for i = 1:1:length(rhoinit2)
        rhoinit2(i) = integral(rho0two,xs(i)-deltax/2,xs(i)+deltax/2);
        rhoinit3(i) = integral(rho0three,xs(i)-deltax/2,xs(i)+deltax/2);
    end
    Hp = @(rho) Hpfunct2(rho,nu,m);
    [entropy2,rho2] = evolve1DEntropy(Hp,@V0,Wconv,xs,deltax,rhoinit2,theta,deltat,k,H);
    [entropy3,rho3] = evolve1DEntropy(Hp,@V0,Wconv,xs,deltax,rhoinit3,theta,deltat,k,H);
    Ts = 0:deltat:Tmax;
end
end

function [entropies,outRhos] = evolve1DEntropy(Hp,V,Wconv,x,deltax,rho0,theta,deltat,k,H)
Potential = V(x)';
 %If x = 0:0.1:1, row1 = 0 - x -> Wconv * rho = desired output 
N = length(rho0);
entropies = zeros(k+1,1);
outRhos = zeros(k+1,N);
outRho = rho0';
outRhos(1,:) = outRho;
entropies(1) = entropy(outRho,deltax,Wconv,H,Potential);
for i = 1:1:k
    outRho=RKStep(Hp,Potential,V,Wconv,x,deltax,outRho,N,theta,deltat);
    outRhos(i+1,:) = outRho;
    entropies(i+1) = entropy(outRho,deltax,Wconv,H,Potential);
    [i/k]
end
end

function [outRho] = RKStep(Hp,Potential,V,Wconv,x,deltax,rho,N,theta,deltat)
    rho1 = rho + deltat*Fluxes1D(Hp,Potential,V,Wconv,x,deltax,rho,N,theta);
    rho2 = 0.75 * rho + 0.25 * rho1 + 0.25 * deltat * Fluxes1D(Hp,Potential,V,Wconv,x,deltax,rho1,N,theta);
    outRho = (rho/3) + (2 * rho2/3) + (2 * deltat * Fluxes1D(Hp,Potential,V,Wconv,x,deltax,rho2,N,theta)/3);
end

function [drhodt] = Fluxes1D(Hp,Potential,V,Wconv,x,deltax,rho,N,theta)
    xij = deltax * Wconv * rho + Hp(rho) + Potential;
    uj = (xij(1:N-1)-xij(2:N))/deltax;
    ujp = max(uj,0);
    ujm = min(uj,0);
    [centeredSlopes,minmodSlopes] = rhox1d(Hp,V,x,deltax,rho,theta);
    rhoEcentered = rho + deltax * centeredSlopes / 2;
    rhoWcentered = rho - deltax * centeredSlopes / 2;
    rhoEminmod = rho + deltax * minmodSlopes / 2;
    rhoWminmod = rho - deltax * minmodSlopes / 2;
    signs = max(sign(rhoEcentered) .* sign(rhoWcentered),0);
    rhoE = rhoEminmod.*(1-signs) + signs .* rhoEcentered;
    rhoW = rhoWminmod.*(1-signs) + signs .* rhoWcentered;
    Fjp = rhoE(1:N-1) .* ujp + rhoW(2:N) .* ujm; 
    Fjplus = [0; Fjp;  0];
    drhodt = (Fjplus(1:N)-Fjplus(2:N+1))/deltax;
end

function [centeredSlopes,minmodSlopes] = rhox1d(Hp,V,x,deltax,rho,theta)
rholeft = ghostCell(Hp,V,x(1),deltax,rho(1));
rhoright = ghostCell(Hp,V,x(length(x)),-deltax,rho(length(rho)));
rhosize = size(rho);
if rhosize(2) == 1
    extendedrho = zeros(rhosize(1)+2,1);
else
    extendedrho = zeros(rhosize(1)+2,rhosize(2)+2);
end
extendedrho(1,:)=rholeft;
extendedrho(2:rhosize(1)+1,:)=rho;
extendedrho(rhosize(1)+2,:)=rhoright;
centeredSlopes = (extendedrho(3:rhosize(1)+2,:) - extendedrho(1:rhosize(1),:))/(2*deltax);
rightSlopes = (extendedrho(3:rhosize(1)+2,:) - extendedrho(2:rhosize(1)+1,:))/deltax;
leftSlopes = (extendedrho(2:rhosize(1)+1,:) - extendedrho(1:rhosize(1),:))/deltax;
minmodSlopes = minmod(theta * rightSlopes,centeredSlopes,theta*leftSlopes);
end

function [rhoOut] = ghostCell(Hp,V,x,deltax,rho)
target = Hp(rho) + V(x) - V(x-deltax);
targetfunc = @(r) Hp(r) - target;
rhoOut = fzero(targetfunc,rho);
end

%Zero Funcs

function y = W0(j)
y = zeros(size(j));
end

function y = H0(j)
y = zeros(size(j));
end

function y = V0(j)
y = zeros(size(j));
end

%Diffusion with Nonlocal attraction
function y = Hpfunct2(x,v,m)
    y = v * x .^ (m-1);
end

%Diffusion with Nonlocal attraction
function y = Hfunct2(x,v,m)
    y = (v/m) * x .^ m;
end

function y = entropy(rho,deltax,Wconv,H,Potential)
    y = deltax*(((deltax/2) *rho' * Wconv * rho)+sum(H(rho))+sum(rho .* Potential));
end

function y = Wfunct2(x)
    y = -max(0,1-abs(x));
end

function y = rho02(x,size)
    y = max(0,sign(size-abs(x)))/size;
end



