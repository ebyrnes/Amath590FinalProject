function [centeredSlopes,minmodSlopes] = rhoy2d(Hp,V,domain,deltay,rho,theta)
ghostfuncbottom = @(rhoval,coords) ghostCell(Hp,V,coords,deltay,rhoval);
ghostfunctop = @(rhoval,coords) ghostCell(Hp,V,coords,-deltay,rhoval);
rhotop = arrayfun(ghostfuncbottom,rho(:,1),domain(:,1));
rholen = length(domain(1,:));
rhobottom = arrayfun(ghostfunctop,rho(:,rholen),domain(:,rholen));
rhosize = size(rho);
extendedrho = zeros(rhosize(1),rhosize(2)+2);
extendedrho(:,1)=rhotop;
extendedrho(:,2:rhosize(1)+1)=rho;
extendedrho(:,rhosize(1)+2)=rhobottom;
centeredSlopes = (extendedrho(:,3:rhosize(1)+2) - extendedrho(:,1:rhosize(1)))/(2*deltay);
rightSlopes = (extendedrho(:,3:rhosize(1)+2) - extendedrho(:,2:rhosize(1)+1))/deltay;
leftSlopes = (extendedrho(:,2:rhosize(1)+1) - extendedrho(:,1:rhosize(1)))/deltay;
minmodSlopes = minmod(theta * rightSlopes,centeredSlopes,theta*leftSlopes);
end

function [rhoOut] = ghostCell(Hp,V,coords,deltay,rho)
target = Hp(rho) + V(coords) - V(coords - 1i* deltay);
targetfunc = @(r) Hp(r) - target;
rhoOut = fzero(targetfunc,rho);
end
