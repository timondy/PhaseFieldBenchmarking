% Code CH2D
%
% Version 1.2 written October 16, 2018 by Brian Wetton, wetton@math.ubc.ca
% Changed initial conditions to the second benchmark problem in 
% "High Accuracy Benchmark Problems for Allen-Cahn and Cahn-Hilliard
%   Dynamics"
% 
% Computes approximations to 2D Cahn-Hilliard equations using the
% methods described in the JCP paper,
%
% "High accuracy solutions to energy gradient flows
%   from material science models"
% by Andrew Christlieb, Jaylan Jones, Keith Promislow, Brian Wetton
%   and Mark Willoughby
%
% Uses subroutines:
%
% uinit:    defines intital data for the computation
% uplot:    outputs a contour plot of solutions
% frhs:     right hand side of the method of lines ODE
% chstep:   calls ch1step, allows generalization to other times stepping
% ch1step:  calls chcg, Newton iterations for implicit substeps
% chcg:     preconditioned conjugate gradient iterations for ch1step
% bigE:     error estimation
% energy:   calculates underlying energy of a solution
%
% Structures param and symbol contain information passed to subroutines
%
%----------
%
% Add timing
%
tic

%
% User defined parameters
%

% set numerical, physical and run parameters
N=192;          % grid resolution 
epsilon=0.1;  % order parameter for the PDE
finalT = 30;   % final time for the comptation 
sigma = 1e-4;   % Local error tolerance 

% Set to 1 to generate a movie
make_movie = 0;
movie_name = 'run.avi';

%----------
% Internal computational parameters

h = 2*pi/N;
param.h=h;
param.N = N;
param.epsilon = epsilon;
param.contourv = linspace(-1.1,1.1,12); % contour lines for the plots 

% maximum total number of CG or fixed point interations per time step
param.maxCG = 500;

% maximum number of Newton steps per time step
param.maxN = 8;

% Minimum and starting time step
mink = 1e-12;
startk = 1e-4;

% Tolerances
param.Ntol = 1e-7;  % Newton iteration residual tolerances 
param.cgtol = 1e-9; % PCG iteration residual tolerances 

% Standard numerical parameters for the adaptive time stepping 
ksafety = 0.8;
kfact = 1.3;
kfact2 = 1/1.3;
Nfact = 0.7;
CGfact = 0.7;

%----------
% Initialize symbols of operators

x = (1:N)*h;
alpha_alias = [0:N/2-1 N/2 -N/2+1:-1];
alpha_alias2 = [0:N/2-1 0 -N/2+1:-1];

xx = zeros(N,N);
yy = zeros(N,N);
alpha = zeros(N,N);
beta = zeros(N,N);
alpha2 = zeros(N,N);
beta2 = zeros(N,N);

for i=1:N
    xx(i,:) = x;
    yy(:,i) = x';
    alpha(i,:) = alpha_alias;
    beta(:,i) = alpha_alias';
    alpha2(i,:) = alpha_alias2;
    beta2(:,i) = alpha_alias2';
end

% x and y coordinates for initial conditions and plotting 
param.xx = xx;
param.yy = yy;

% initial conditions
u0 = uinit(param);

% operator symbols
lap = - (alpha.^2+beta.^2);
lap2 = lap;
% lap2 is the symbol of the Laplacian on massless solutions
% the value below is arbitrary but avoids 0/0 in the H_{-1} PCG inner products
lap2(1,1) = 1;
bih = lap.*lap;

symbol.lap = lap;
symbol.lap2 = lap2;
symbol.biharmonic = bih;
symbol.derivx = sqrt(-1)*alpha2;
symbol.derivy = sqrt(-1)*beta2;

if make_movie == 1
    video = VideoWriter(movie_name);
    open(video);
end

%----------
% plot initial condition
uplot(u0,param,1,'Initial Conditions');

%----------
% initialize Energy and iteration counters
% and solution values at (pi/2,pi/2) and (3pi/2,3pi/2)

fu0 = frhs(u0,param,symbol);
E0 = energy(u0,param,symbol);
E = E0;
No4 = N/4;

U1 = u0(No4,No4);
U2 = u0(3*No4,3*No4);
% Flag for zero crossings for the benchmark times 
U1flag = 0;
U2flag = 0;

Nit = [];
CGit = [];

it_tot = 0;
cgit_tot = 0;

tt = 0;
t = tt;
kk = startk;
k = [];
Nt = 1;

%----------
% Main loop through the adaptive time stepping

epic_fail = 0; % if things really cannot be continued 
disp('Begin computation, reporting every 100 time steps')

while tt < finalT && epic_fail==0
    % one DIRK2 time step
    [u1 fail nloc cgloc] = chstep(kk,param,u0,fu0,symbol);
    if fail==0
        fu1 = frhs(u1,param,symbol);
        E1 = energy(u1,param,symbol);
        if E1 > E0*(1+sigma)
            disp('Energy increase, failing time step')
            fail = 1;
        end
    end
    
    % check accuracy if time step computation was successful 
    accfail = 0;
    if fail == 0
        Xi = bigE(u0,u1,fu0,fu1,kk, param, symbol);
        if  Xi > sigma
            % disp('Accuracy lost, failing time step')
            accfail =1;
        end
    end
    
    % check if time step failed or was innacurate 
    if fail ==1 || accfail ==1
        if fail ==1 || Xi/sigma > 2
            kk = kk*kfact2;
        else
            kk = kk*ksafety*sqrt(sigma/Xi);
        end
        
        if kk<mink
            if fail ==1
                epic_fail =1;
                disp('Failed time step fell below minimum -> epic fail')
            else
                kk = mink;
                disp('Minimum time step used')
            end
        end
    else   
    % accept time step
    % bookkeep step and iteration counts
        
        Nit = [Nit nloc];
        CGit = [CGit cgloc];
        k = [k kk];
        it_tot = it_tot+nloc;
        cgit_tot = cgit_tot + cgloc;
        t = [t tt+kk];
        tt = tt+kk;
        Nt = Nt+1;
        E = [E E1];
        
        % Check for first benchmark time
        % Use linear interpolation between the time steps for the value 
        U1 = [U1 u0(No4,No4)];
        if U1flag == 0 && u0(No4,No4) < 0 
            U1flag = 1;
            t2 = t(end);
            v2 = U1(end);
            t1 = t(end-1);
            v1 = U1(end-1);
            ttrans = t1+v1/(v1-v2)*(t2-t1);
            U1ttrans = ttrans;
            fprintf('U1 changed sign at time %d \n',ttrans)
        end
        
        % Check for second benchnark time 
        U2 = [U2 u0(3*No4,3*No4)];
        if U2flag == 0 && u0(3*No4,3*No4) < 0 
            U2flag = 1;
            t2 = t(end);
            v2 = U2(end);
            t1 = t(end-1);
            v1 = U2(end-1);
            ttrans = t1+v1/(v1-v2)*(t2-t1);
            U2ttrans = ttrans;
            fprintf('U2 changed sign at time %d \n',ttrans)
        end
        
        % period plot every 100 time steps 
        if mod(Nt,100) == 0
            fprintf('Advanced to time %d with time step %d \n',tt,kk)
            if make_movie == 0
                uplot(u0,param,7,'Periodic plot');
            end
        end
        u0 = u1;
        fu0 = fu1;
        E0 = E1;
        
        % plot time steps to movie
        if make_movie == 1
            uplot(u0,param,8, ...
                sprintf('ln k = %8.2f,   t = %10.3e',log(kk),tt)); 
            writeVideo(video,getframe(gcf));
        end
        
        % estimate appropriate next time step
        kknew = kk*ksafety*sqrt(sigma/Xi);
        kknew = min(kknew, kk*kfact);
        
        % don't increase time step if iteration counts are too close to tolerance
        if nloc/param.maxN > Nfact || cgloc/param.maxCG > CGfact
            kknew = min(kknew, kk);
        end
        
        kk = kknew;
        
    end
end

if epic_fail == 1
    disp('Time step reduced below minimum')
    disp('*** Computation failed ***')
else
    % successful run, output iteration counts and graphs
    fprintf('Total CG iterations %d \n',cgit_tot)
    fprintf('Total time steps %d \n \n', Nt)
    
    uplot(u0,param,2,'Final Solution');
    
    figure(3)
    plot(t,E)
    title('Energy')
    
    figure(4)
    plot(t(1:Nt-1), k)
    title('Time steps k')
    
    figure(5)
    plot(t(1:Nt-1), Nit)
    title('Newton iterations')
    
    figure(6)
    plot(t(1:Nt-1), CGit)
    title('Conjugate gradient iterations')
    
    figure(7) 
    plot(t,U1,t,U2)
    legend('(\pi/2,\pi/2)', '(3\pi/2,3\pi/2)')
    title('u values')
    
    if make_movie ==1
       close(video);
    end
end

% Save results
save('run.mat','U1','U2','E','t','U1ttrans','U2ttrans')

% Output computational time 
fprintf('\n Computational time: %d \n',toc)



