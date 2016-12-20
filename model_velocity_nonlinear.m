function [unl, xC0, yC0, angC0] = ...
    model_velocity_nonlinear(xC,yC,kC,Qtot,B,H,U,Cf,slope,nonlinear,plt)
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
iR0 = -kC;
angC0 = angle((yC1-1i*xC1));
o = ones(size(sC0));


% [unl] = meander_model_AR_copy(sC0,iR0,0*o,Qtot/B*o,H*o,B*o,Cf,slope,-H*o,0*o,0.01,nonlinear,0,1,1,0.8,0.25);
%  for k=1:10
%      [unl] = meander_model_AR_copy(unl.s,iR0,unl.AR,unl.q,unl.h,unl.W,Cf,slope,-H*o,unl.asR,0.01,nonlinear,0,1,1,0.8,0.25);
%  end

[unl] = meander_model_copy(sC0,iR0,0,Qtot/B*o,H*o,B*o,Cf,slope,-H*o,0*o,0.01,nonlinear,0);
for k = 1:10
    [unl] = meander_model_copy(unl.s,iR0,unl.AR,unl.q,unl.h,unl.W,Cf,slope,-H*o,unl.asR,0.01,nonlinear,0);
end

% [ul] = meander_model_copy(sC0,iR0,4*o,Qtot/B*o,H*o,B*o,Cf,slope,-H*o,0*o,0.01,0,0);
% for k=1:10
%     [ul] = meander_model_copy(ul.s,iR0,ul.AR,ul.q,ul.h,ul.W,Cf,slope,-H*o,ul.asR,0.01,0,0);
% end




end


