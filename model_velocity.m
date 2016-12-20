function [ub] = model_velocity(xi,yi,s,k,W,H,U,Cf,psi)

    rho = 998.1;    % density of water, kg/m^3
    g = 9.81;       % gravity, m/s^2
    B = W./2;
    aspect_ratio = W./H;
    alphaa = 10;
    betaao = B/H;           % half-width to depth ratio for straight channel characterized by valley slope
    F2 = (U/sqrt(g*H))^2;   % initial (Froude #)^2 for straight channel characterized by valley slope
    appx_conv_int = 0;      % approximate convolution integral by truncating its computation after some threshold distance (currently hard-coded to be 100B)
    
    %%% scour factor
    %Constantine et al 2009, assuming max depth at bank toe
    % Ikeda 1981: depth at toe relative to mean depth is prop to curvature.
    % A is the proportionality coeff. A = rel. depth ./ curvature
    A = 3.79; % based on linear fitting field data from sites, using Ikeda 1981 fn
  
    A = 10; %approximate maximum from field obs
    chi1 = 0.077./sqrt(Cf); % J&P 1988
    chi = chi1 - 1/3;
    As = 181*(H/B)^2/chi1*(2*chi^2 + 4/5*chi + 1/15);
    alphaa = A + As - 1; % Constantine et al 2009
    %%%
    ub = flowfield_Schwenk(B,k',H,s',alphaa,F2,Cf,appx_conv_int);
    %ub = flowfield_Schwenk_modified(B,k',H,s',alphaa,F2,Cf,appx_conv_int,psi);
    % ub: dimensionless velocity perturbation (fraction of avg velocity)

end

