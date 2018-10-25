function u0 = uinit(param)
% initial conditions for the run

xx = param.xx;
yy = param.yy;
epsilon = param.epsilon;
N = param.N;

% Bechmark problem #2 starts with seven smoothed disks of phase +1

u0 = -ones(N,N);

Ncirc = 7;
xcirc = [1/4 1/8 1/4 1/2 3/4 1/2 3/4]*2*pi;
ycirc = [1/4 3/8 5/8 1/8 1/8 1/2 3/4]*2*pi;
rcirc = [1/10 1/15 1/15 1/20 1/20 1/8 1/8]*2*pi;

for i=1:N
    for j=1:N
        for c=1:Ncirc
            r = sqrt((xx(i,j)-xcirc(c))^2 + (yy(i,j)-ycirc(c))^2);
            if r - rcirc(c) < 0  
                u = -1+2*exp(-epsilon^2/(r-rcirc(c))^2);
                u0(i,j) = u;
            end
        end
    end
end


