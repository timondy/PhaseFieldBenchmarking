function [xk,fail,iter] = chcg(k,param,u,r0,alpha,symbol)
% Preconditioned CG gradient iterations for implicit steps 
% Standard except for the use of the H^{-1} inner product 

fail =0;
N = param.N;
lap = symbol.lap;
lap2 = symbol.lap2;
epsilon = param.epsilon;
bih = symbol.biharmonic;
maxcg = param.maxCG;
cgtol = param.cgtol;

q = 1+alpha*k*bih*epsilon^2-2*alpha*k*lap;

xk = zeros(N,N);
rk = r0;
zk = real(ifft2(fft2(rk)./q));
pk = zk;
linvzk = -real(ifft2(fft2(zk)./lap2));
ipnew = sum(sum(linvzk.*rk));

resid = 1;
iter = 0;

while resid>cgtol && iter <maxcg
    ipold = ipnew;
    Apk1 = -alpha*k*real(ifft2(fft2((3*u.^2-1).*pk).*lap));
    Apk2 = alpha*k*epsilon^2*real(ifft2(fft2(pk).*bih));
    Apk = pk+Apk1+Apk2;
    % Laplacian inverse needed here for H_{-1} inner products
    linvpk = -real(ifft2(fft2(pk)./lap2));
    cgalpha = ipold/sum(sum(linvpk.*Apk));
    xk = xk+cgalpha*pk;
    rk = rk - cgalpha*Apk;
    zk = real(ifft2(fft2(rk)./q));
    linvzk = -real(ifft2(fft2(zk)./lap2));
    ipnew = sum(sum(rk.*linvzk));
    beta = ipnew/ipold;
    pk = zk+beta*pk;
    resid = max(max(abs(rk)));
    iter = iter+1;
end

if resid>cgtol
    disp('Too many CG steps')
    fail =1;
end