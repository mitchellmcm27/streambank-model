function [unl, xC0, yC0, angC0] = ...
    model_velocity_nonlinear_swe(xC,yC,kC,Q,B,H,sC,Cf,slope,nonlinear,AR0,bedtype)
xC=xC(:)';
yC=yC(:)';

%% Savitsky-Golay filter settings
degree = 5; % 5
halfwin = 7; % 7
intnum = 2; % 2

%% curvature
sC=[0,cumsum(sqrt(diff(xC).^2+diff(yC).^2))];
[xC0,xC1,xC2,xm0] = sgolayirreg(sC,xC,degree,halfwin);
[yC0,yC1,yC2,ym0] = sgolayirreg(sC,yC,degree,halfwin);

sC0=sC;
zt = (xC1.^2 + yC1.^2);
iR0 = sign(yC1.*xC2 - xC1.*yC2).*sqrt(zt.*(xC2.^2 + yC2.^2)-((xC1.*xC2+yC1.*yC2).^2))./(zt.^(3/2));
angC0 = angle((yC1-1i*xC1));
o = ones(size(sC0));
zb = 0*o;

[unl] = meander_model_copy(sC,iR0,AR0,Q/B*o,H*o,B*o,Cf,slope,-H*o,o*0,0.01,nonlinear,0,bedtype);
for k = 1:3
    [unl] = meander_model_copy(unl.s,iR0,AR0,unl.q,unl.h,unl.W,Cf,slope,-H*o,unl.asR,0.01,nonlinear,0,bedtype);
end
end


