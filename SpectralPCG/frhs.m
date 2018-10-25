function f = frhs(u, param, symbol)
% computes the RHS of the time dependent problem 

f = -param.epsilon^2*real(ifft2(fft2(u).*symbol.biharmonic)) + ...
    real(ifft2(fft2(u.^3-u).*symbol.lap));