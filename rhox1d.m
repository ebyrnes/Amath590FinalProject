

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
