function [unew,fail,nloc,cgloc] = chstep(k, param, u0, fold, symbol)
% Wrapper for ch1step that allows for easy implementation of substeps in
% higher order methods

% Here, DIRK2 is implemented 

alpha = (2-sqrt(2))/2;

[fstar, fail, nloc1, cgloc1] = ch1step(k, param, u0, alpha, fold, symbol);
ustar = u0 + k*(1-alpha)*fstar;

if fail == 0
    [fnew, fail, nloc2, cgloc2] = ch1step(k, param, ustar, alpha, fstar, symbol);
    unew = ustar + k*alpha*fnew;
else
    nloc2=0;
    cgloc2=0;
    unew = u0;
end

nloc = nloc1+nloc2;
cgloc = cgloc1+cgloc2;