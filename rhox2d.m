function [centeredSlopes,minmodSlopes] = rhox2d(Hp,V,domain,deltax,rho,theta)
ghostfuncleft = @(rhoval,coords) ghostCell(Hp,V,coords,deltax,rhoval);
ghostfuncright = @(rhoval,coords) ghostCell(Hp,V,coords,-deltax,rhoval);
rholeft = arrayfun(ghostfuncleft,rho(1,:),domain(1,:));
rholen = length(domain(:,1));
rhoright = arrayfun(ghostfuncright,rho(rholen,:),domain(rholen,:));
rhosize = size(rho);
extendedrho = zeros(rhosize(1)+2,rhosize(2));
extendedrho(1,:)=rholeft;
extendedrho(2:rhosize(1)+1,:)=rho;
extendedrho(rhosize(1)+2,:)=rhoright;
centeredSlopes = (extendedrho(3:rhosize(1)+2,:) - extendedrho(1:rhosize(1),:))/(2*deltax);
rightSlopes = (extendedrho(3:rhosize(1)+2,:) - extendedrho(2:rhosize(1)+1,:))/deltax;
leftSlopes = (extendedrho(2:rhosize(1)+1,:) - extendedrho(1:rhosize(1),:))/deltax;
minmodSlopes = minmod(theta * rightSlopes,centeredSlopes,theta*leftSlopes);
end

function [rhoOut] = ghostCell(Hp,V,coords,deltax,rho)
target = Hp(rho) + V(coords) - V(coords - deltax);
targetfunc = @(r) Hp(r) - target;
rhoOut = fzero(targetfunc,rho);
end
