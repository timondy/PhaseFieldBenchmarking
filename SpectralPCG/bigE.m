function E = bigE(u0,u1,f0,f1,k,param,symbol)
% estimates the local error of a times step using a high order predictor 

uhalf = (u0+u1)/2+k*(f0-f1)/8;
fhalf = frhs(uhalf,param,symbol);
uaccurate = u0+k*(f0 + 4*fhalf + f1)/6;
E = max(max(abs(u1-uaccurate)));
