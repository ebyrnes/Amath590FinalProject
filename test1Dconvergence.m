function [outData] = test1Dconvergence(h)
if h == 1
    X = 2;
    Nmin = 3;
    Nmax = 9;
    outData = zeros(3,Nmax+1-Nmin);
    for j = Nmin:1:Nmax
        deltax = 2^(2-j);
        deltat = deltax/2;
        k = ceil(750/deltat);
        xs = -X:deltax:X;
        rhoinit = zeros(size(xs));
        rhoinf = zeros(size(xs))';
        for i = 1:1:length(rhoinit)
            rhoinit(i) = integral(@rho01,xs(i)-deltax/2,xs(i)+deltax/2);
            rhoinf(i) = integral(@rhoinf1,xs(i)-deltax/2,xs(i)+deltax/2);
        end
        rhoinit = rhoinit * sum(rhoinf)/sum(rhoinit);
        theta=2;
        Wconv = Wfunct1(xs'-xs,deltax);
        outRho = evolve1DSS(@H0,@V0,Wconv,xs,deltax,rhoinit,theta,deltat,k,j);
        L1Norm = deltax*sum(abs(outRho-rhoinf));
        LinfNorm = max(abs(outRho-rhoinf));
        outData(1,1+j-Nmin) = deltax;
        outData(2,1+j-Nmin) = L1Norm;
        outData(3,1+j-Nmin) = LinfNorm;
    end
end
end

function [outRho] = evolve1DSS(Hp,V,Wconv,x,deltax,rho0,theta,deltat,k,j)
Potential = V(x)';
 %If x = 0:0.1:1, row1 = 0 - x -> Wconv * rho = desired output 
N = length(rho0);
outRho = rho0';
for i = 1:1:k
    outRho=RKStep(Hp,Potential,V,Wconv,x,deltax,outRho,N,theta,deltat);
    [i/k, j]
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

%Attraction-repulsion kernel

function y = rho01(x)
y = exp(-x.^2/2)/sqrt(2*pi);
end

function y = rhoinf1(x)
y = real(sqrt(2-x.^2)/pi);
end

function y = Wfunct1(x,deltax)
abX = abs(x);
y = ((abX .^2 / 2) - log(abX));
y(1:size(y,2)+1:end)=(1+(deltax^2/24)-log(deltax/2));
end


