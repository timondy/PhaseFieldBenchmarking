function T = computeT(fname)

dat = load(fname);

t = dat(:,1);
u = dat(:,2);

i0=find(u(2:end).*u(1:end-1)<=0);
if (i0)
    i1=i0+1;
    T=t(i0) - (t(i1)-t(i0)) * u(i0)/(u(i1)-u(i0));
else
    T=NaN;
end

end