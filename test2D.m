function outRhos = test2D(Tmax)
X = 4;
Y = 4;
deltax = 0.1;
deltay = 0.1;
xs = -X:deltax:X;
ys = -Y:deltay:Y;
Nx = length(xs);
Ny = length(ys);
coords = xs' + 1i * ys;
xconv = xs'-xs;
yconv = ys'-ys;
domainconv = zeros(length(xs),length(xs),length(ys),length(ys));
for i = 1:1:length(ys)
    for j = 1:1:length(ys)
        domainconv(:,:,i,j) = xconv + 1i * yconv(i,j);
    end
end
Wconv = W1(domainconv);
deltat = 0.001;
nsteps = Tmax / deltat;
v = 0.1;
theta = 2;
Hp = @(x) Hpfunct2(x,v,3);
domainlim = 3;
rho0f = @(x,y) rho0func2(x,y,domainlim,domainlim);
rho0 = zeros(size(coords));
for i = 1:1:length(xs)
    for j = 1:1:length(ys)
        rho0(i,j)=integral2(rho0f,xs(i)-(deltax/2),xs(i)+(deltax/2),ys(j)-(deltay/2),ys(j)+(deltay/2))/(deltax * deltay);
    end
end
outRhos = evolve2D(Hp,@V0,Wconv,coords,deltax,deltay,rho0,theta,deltat,nsteps,Nx,Ny);
end

function [outRhos] = evolve2D(Hp,V,Wconv,coords,deltax,deltay,rho0,theta,deltat,nsteps,Nx,Ny)
Potential = V(coords);
sz = size(rho0);
outRhos = zeros(nsteps+1,sz(1),sz(2));
outRhos(1,:,:) = rho0;
for i = 1:1:nsteps
    i * deltat
    rho = reshape(outRhos(i,:,:),Nx,Ny);
    rho1 = rho + deltat*Fluxes2D(Hp,Potential,V,Wconv,coords,deltax,deltay,rho,Nx,Ny,theta);
    rho2 = 0.75 * rho + 0.25 * rho1 + 0.25 * deltat * Fluxes2D(Hp,Potential,V,Wconv,coords,deltax,deltay,rho1,Nx,Ny,theta);
    outRhos(i+1,:,:)= (rho/3) + (2 * rho2/3) + (2 * deltat * Fluxes2D(Hp,Potential,V,Wconv,coords,deltax,deltay,rho2,Nx,Ny,theta)/3);
end

end

function [drhodt] = Fluxes2D(Hp,Potential,V,Wconv,coords,deltax,deltay,rho,Nx,Ny,theta)
    convTerm = deltax * deltay * tensorprod(Wconv,rho,[2,4],[1,2]);
    xij = convTerm + Hp(rho) + Potential;
    uj = (xij(1:Nx-1,:)-xij(2:Nx,:))/deltax;
    ujp = max(uj,0);
    ujm = min(uj,0);
    vj = (xij(:,1:Ny-1)-xij(:,2:Ny))/deltay;
    vjp = max(vj,0);
    vjm = min(vj,0);
    [centeredSlopesx,minmodSlopesx] = rhox2d(Hp,V,coords,deltax,rho,theta);
    rhoEcentered = rho + deltax * centeredSlopesx / 2;
    rhoWcentered = rho - deltax * centeredSlopesx / 2;
    rhoEminmod = rho + deltax * minmodSlopesx / 2;
    rhoWminmod = rho - deltax * minmodSlopesx / 2;
    signsEW = max(sign(rhoEcentered) .* sign(rhoWcentered),0);
    rhoE = rhoEminmod.*(1-signsEW) + signsEW .* rhoEcentered;
    rhoW = rhoWminmod.*(1-signsEW) + signsEW .* rhoWcentered;

    Fjp = rhoE(1:Nx-1,:) .* ujp + rhoW(2:Nx,:) .* ujm; 
    Fjplus = [zeros(1,length(Fjp(1,:))); Fjp; zeros(1,length(Fjp(1,:)))];
    drhodt1 = (Fjplus(1:Nx,:)-Fjplus(2:Nx+1,:))/deltax;

    [centeredSlopesy,minmodSlopesy] = rhox2d(Hp,V,coords,deltay,rho,theta);
    rhoNcentered = rho + deltay * centeredSlopesy / 2;
    rhoScentered = rho - deltay * centeredSlopesy / 2;
    rhoNminmod = rho + deltay * minmodSlopesy / 2;
    rhoSminmod = rho - deltay * minmodSlopesy / 2;
    signsNS = max(sign(rhoNcentered) .* sign(rhoScentered),0);
    rhoN = rhoNminmod.*(1-signsNS) + signsNS .* rhoNcentered;
    rhoS = rhoSminmod.*(1-signsNS) + signsNS .* rhoScentered;

    Fkp = rhoN(:,1:Ny-1) .* vjp + rhoS(:,2:Ny) .* vjm; 
    Fkplus = [zeros(length(Fkp(:,1)),1) Fkp zeros(length(Fkp(:,1)),1)];
    drhodt2 = (Fkplus(:,1:Ny)-Fkplus(:,2:Ny+1))/deltay;

    drhodt = drhodt1 + drhodt2;
end



function out = W1(z)
out = -exp(-abs(z).^2)/pi;
end

function y = V0(j)
y = zeros(size(j));
end

function y = Hpfunct2(x,v,m)
    y = v * x .^ (m-1);
end

function out = rho0func(x,xlim,ylim)
    if (abs(real(x)) <= xlim) && (abs(imag(x)) <= ylim)
        out = 1/4;
    else
        out = 0;
    end
end

function out = rho0func2(x,y,xlim,ylim)
    out = min(max(0,sign(xlim-abs(x))),max(0,sign(ylim-abs(y))))/4;
end