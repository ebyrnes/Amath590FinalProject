function [out] = minmod(x,y,z)
% Computes minmod in a vectorized way.
maxes = max(max(x,y),z);
mins = min(min(x,y),z);
maxsigns = sign(maxes);
minsigns = sign(mins);
signvals = (maxsigns + minsigns)/2;
posPart = (1 + signvals) .* mins / 2;
negPart = abs((signvals - 1)) .* maxes / 2;
out = abs(signvals) .* (posPart + negPart);
end