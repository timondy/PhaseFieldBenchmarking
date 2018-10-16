function [fnew,fail,nloc,cgloc] = ch1step(k, param, ustar, alpha, ...
    fold, symbol)
% Newton method iterations for implicit substeps of fully discrete CH
% This is one of the two substeps for the DIRK2 time stepping, with
%   parameter alpha 

% Report failure is iteration count is exceeded...
%   or Newton iterations go "wild"  
fail = 0;
maxN = param.maxN;
maxcg = param.maxCG;
Ntol = param.Ntol;

cgloc = 0;
nloc = 0;
resid = 1;

fnew = fold;

% Newton method to reduce the residual to the nonlinear equations 
rhs = fnew - frhs(ustar+k*alpha*fnew,param,symbol);

while resid>Ntol && nloc < maxN && fail == 0
    [delta fail cgit] = chcg(k,param,ustar+k*alpha*fnew,rhs,alpha,symbol);
    cgloc = cgloc +cgit;
    nloc = nloc +1;
    fnew = fnew - delta;
    if max(max(abs(k*fnew))) > 5
        fail = 1;
        disp('wild Newton Step')
    else
        rhs = fnew - frhs(ustar+k*alpha*fnew,param,symbol);
        resid = max(max(abs(rhs)));
    end
end

if resid>Ntol
    disp('maximum iterations exceeded')
    fail =1;
end

if cgloc>maxcg
    disp('maximum CG iterations exceeded')
    fail = 1;
end

