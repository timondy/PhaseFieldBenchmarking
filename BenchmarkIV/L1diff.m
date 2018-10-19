% 
% Compute the L1 error in [ln(E)](theta) where theta is ln(t) between 
% two Benchmark IV runs.
% 

% Substitute your data for run8
load('run8.mat')
lt1 = log(t(2:end));
lE1 = log(E(2:end));

% The most accurate result we have 
load('run9.mat')
lt2 = log(t(2:end));
lE2 = log(E(2:end));

theta = linspace(-5,7,1000); % Use for benchmark D1
% theta = linspace(-5,2,1000); % Use for benchmark D2

% Linearly interpolate values to a uniform grid
lE1interp = interp1(lt1,lE1,theta);
lE2interp = interp1(lt2,lE2,theta);

diff = trapz(theta, abs(lE1interp-lE2interp));
fprintf('D value: %d \n', diff)