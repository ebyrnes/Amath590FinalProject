function [outRhos] = evolve1D(Hp,V,Wconv,x,deltax,rho0,theta,deltat,nsteps)
Potential = V(x)';
 %If x = 0:0.1:1, row1 = 0 - x -> Wconv * rho = desired output 
N = length(rho0);
outRhos = zeros(nsteps+1,N);
outRhos(1,:) = rho0;
for i = 1:1:nsteps
    outRhos(i+1,:)=RKStep(Hp,Potential,V,Wconv,x,deltax,outRhos(i,:)',N,theta,deltat);
end
end

function [outRho] = RKStep(Hp,Potential,V,Wconv,x,deltax,rho,N,theta,deltat)
    rho1 = rho + deltat*Fluxes1D(Hp,Potential,V,Wconv,x,deltax,rho,N,theta);
    rho2 = 0.75 * rho + 0.25 * rho1 + 0.25 * deltat * Fluxes1D(Hp,Potential,V,Wconv,x,deltax,rho1,N,theta);
    outRho = (rho/3) + (2 * rho2/3) + (2 * deltat * Fluxes1D(Hp,Potential,V,Wconv,x,deltax,rho2,N,theta)/3);
end