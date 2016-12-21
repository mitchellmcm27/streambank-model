function [out_x, out_y, new_s, k, t, n] = ParaSpline (dpx, dpy, varargin)
%ParaSpline takes "digitzed points" in x and y (dpx, dy),
% interpolates parametrically between them in the order given,
% creates a piecewise polynomial representation F(X(t),Y(t)),
% calculates curvature using first and 2nd derivatives of X and Y,
% resamples the polynomials given input resolution (default 1),
% outputs x,y,t, and curvature at each resampled point.
%
% Usage:
%   ParaSpline(dpx,dpy)
%   ParaSpline(dpx,dpy,output_res)

switch nargin
    case 2
        resampling_interval = 1; % same units as X and Y (meters)
    case 3
        resampling_interval = varargin{1}; %
    otherwise
        error('Invalid inputs. Needed 2 vectors, X and Y, and optionally output resampling interval.');
end

x = dpx(:);
y = dpy(:);
xy = [x y];
xy = unique(xy,'rows','stable');
x=xy(:,1);
y=xy(:,2);
n = numel(x);
m = numel(y);

assert(m==n, 'X and Y vectors not the same length')

s = [0; cumsum(hypot(diff(x),diff(y)))];

pp_x = interp1(s, x, 'spline', 'pp');
pp_y = interp1(s, y, 'spline', 'pp');

dx_pp = fnder(pp_x);
ddx_pp = fnder(pp_x,2);
dy_pp = fnder(pp_y);
ddy_pp = fnder(pp_y,2);

%new_s = [0:resampling_interval:max(s) s']; % resample but keep original points
new_s = [0:resampling_interval:max(s)]; % resample and do not keep original points
new_s = unique(sort(new_s, 'ascend'));
out_x = ppval(pp_x, new_s); %row vector
out_y = ppval(pp_y, new_s);

dx = ppval(dx_pp, new_s);
dy = ppval(dy_pp, new_s);
ddx = ppval(ddx_pp, new_s);
ddy = ppval(ddy_pp, new_s);

k = (dx.*ddy - dy.*ddx)./...
    power(dx.*dx + dy.*dy, 3/2);
t = atan2(dy,dx);
n = atan2(-dx,dy);

end

