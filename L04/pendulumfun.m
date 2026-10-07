function [fiv] = pendulumfun(t,fi,ell)
%PENDULUMFUN Differentia
%   fi'' =  -7*g/(3*1.76*ell)*cos fi 
% 
n=4;
g = 9.81;
fiv = zeros(size(fi));
fiv(1) = fi(2);
fiv(2) = -g/ell*(2*n-1)/1.88*sin(fiv(1));
end

