function uplot(u, param, frame, ftitle)
% plots the current solution

figure(frame)
[cc hh] = contourf(param.xx,param.yy,u,param.contourv);
colorbar;
title(ftitle)
getframe;