function E = energy(u, param, symbol)
% computes the energy of the current solution 

h = param.h;
uhat = fft2(u);
upx = real(ifft2(symbol.derivx.*uhat));
upy = real(ifft2(symbol.derivy.*uhat));
E = h^2*sum(sum(param.epsilon^2*(upx.^2+upy.^2)/2+(u.^4-2*u.^2+1)/4));