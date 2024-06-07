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
