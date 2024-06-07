function outRhos = test1D(h)

if h ==1
    disp("Attractive-Repulsive kernel")
    X = 2;
    T = 2;
    m = 50;
    k = 1000;
    deltax = 2*X/(2*m+1);
    deltat = T/k;
    xs = linspace(-X,X,((2*m)+1));
    rhoinit = zeros(size(xs));
    rhoinf = zeros(size(xs))';
    for i = 1:1:length(rhoinit)
        rhoinit(i) = integral(@rho01,xs(i)-deltax/2,xs(i)+deltax/2);
        rhoinf(i) = integral(@rhoinf1,xs(i)-deltax/2,xs(i)+deltax/2);
    end
    rhoinit = rhoinit * sum(rhoinf)/sum(rhoinit);
    theta=2;
    Wconv = Wfunct1(xs'-xs,deltax);
    outRhos = evolve1D(@H0,@V0,Wconv,xs,deltax,rhoinit,theta,deltat,k);

elseif h ==2
    disp("Nonlinear Diffusion with nonlocal attraction kernel")
    X = 6;
    T = 100;
    m = 300;
    deltat = 0.001;
    k = T / deltat;
    v = 1.48;
    M = 3;
    sigma = 1;
    deltax = 2*X/(2*m+1);
    xs = linspace(-X,X,((2*m)+1));
    rhoinitfunc = @(x) rho02(xs,M);
    rhoinit = zeros(size(xs));
    for i = 1:1:length(rhoinit)
        rhoinit(i) = integral(rhoinitfunc,xs(i)-deltax/2,xs(i)+deltax/2);
    end
    Hp = @(x) Hpfunct2(x,v,3);
    theta = 2;
    Wconv = Wfunct2(xs'-xs,sigma);
    outRhos = evolve1D(Hp,@V0,Wconv,xs,deltax,rhoinit,theta,deltat,k);
end

%Setting various functions to zero

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
y = (abX .^2 / 2) - log(abX);
y(1:size(y,2)+1:end)=1+(deltax^2/24)-log(deltax/2);
end

%Diffusion with Nonlocal attraction
function y = Hpfunct2(x,v,m)
    y = v * x .^ (m-1);
end

function y = Wfunct2(x,sigma)
    y = -exp((-abs(x).^2)/(2*sigma)) / sqrt(2 * pi * sigma);
end

function y = rho02(x,shift)
    y = (exp(-0.5 * (x - shift).^2)+exp(-0.5 * (x + shift).^2))/sqrt(8*pi);
end


