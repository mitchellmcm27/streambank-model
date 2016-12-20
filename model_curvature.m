function [xi,yi, k, n, s, bank_points] = model_curvature(dpx,dpy,B,site_x,site_y)

dpx = dpx(~isnan(dpx));
dpy = dpy(~isnan(dpy));

resampling_interval = B/10;
[xi, yi, s, k, t, n] = ParaSpline(dpx, dpy, resampling_interval);

% site_cl_idx = dsearchn([xi; yi]', [site_x site_y]);
% new_x = xi(site_cl_idx);
% new_y = yi(site_cl_idx);
%
% point1 = dsearchn([dpx; dpy]',[new_x new_y]);
% point2 = dsearchn([dpx([point1-1 point1+1]); dpy([point1-1 point1+1])]',[new_x new_y]);
% if point2 == 1
%     point2 = point1 - 1;
% elseif point2 == 2
%     point2 = point1 + 1;
% end
% assert(abs(point1-point2)==1);
% 
% if point1>point2
%     dpx = [dpx(1:point2) new_x dpx(point1:end)];
%     dpy = [dpy(1:point2) new_y dpy(point1:end)];
% else
%     dpx = [dpx(1:point1) new_x dpx(point2:end)];
%     dpy = [dpy(1:point1) new_y dpy(point2:end)];
% end
% dp_idx = find(ismember(xi,dpx) & ismember(yi,dpy));

bank_points = [xi+B*sign(k).*cos(n);
    yi+B*sign(k).*sin(n)]';

end

