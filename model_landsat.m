function [vals,dates] = model_landsat(point, pathrow, dates, type)
%% Create a timeseries of Landsat 8 pixel(s).
    function ts = average_monthly_ts(t,x)
        %% Given a daily timeseries, create an average monthly one
        [Y,M] = datevec(t);
        [YM,~,c] = unique([Y' M'],'rows');
        YMval = [YM accumarray(c, x, [], @(x)mean(x))];
        
        a = find(YM(:,1)==year1 & YM(:,2)==month1, 1, 'first');
        b = find(YM(:,1)==year2 & YM(:,2)==month2, 1, 'last');
        if isempty(b)
            b = size(YMval,1);
        end
        val = YMval(a:b,3);
        mos = datetime(year1,month1,15) + calendarDuration(0,find(val)-1,0);
        ts ={mos val};
    end

    function [result, conf] = is_cloudy(evifile,y,x)
        mname = strcat(evifile(1:21),'_cfmask.tif');
        cname = strcat(evifile(1:21),'_cfmask_conf.tif');
        cfmask = fullfile(lsdir,'sr',mname);
        cfmask_conf = fullfile(lsdir,'sr',cname);
        result = imread(cfmask, 'PixelRegion',{[floor(y) ceil(y)], [floor(x) ceil(x)]});
        result = double(result);
        conf = imread(cfmask_conf, 'PixelRegion',{[floor(y) ceil(y)], [floor(x) ceil(x)]});
    end


lsdir = fullfile('data');
interpmethod = 'sin4';
easting = point(1);
northing = point(2);
date1 = dates(1);
date2 = dates(2);
[year1, month1, ~] = datevec(date1);
[year2, month2, ~] = datevec(date2);

d1 = day(date1,'dayofyear');
d2 = day(date2,'dayofyear');

% Set up filenames and z-scales.
switch type
    case 'EVI'
        folder = 'sr';
        fbase = 'LC80';
        fend = 'LGN00_sr_evi';
        ftype = 'tif';
        scale = 0.0001;
        spike_thresh = 0.07/30;
        lims = [-0.1 0.35];
        
    case 'NDVI'
        folder = 'sr';
        fbase = 'LC80';
        fend = 'LGN00_sr_ndvi';
        ftype = 'tif';
        scale = 0.0001;
        spike_thresh = 0.14/30; %slope
        lims = [-0.1 1];
    case 'FC'
        prstring = num2str(pathrow);
        scene = strcat(fullfile(lsdir,'forest',filesep),...
            'p0',prstring(1:2),'r',prstring(3:5),'_TC_2005.tif');
        info = geotiffinfo(scene);
        refmat = info.RefMatrix;
        %pixel_offset = -15;
        pixel_offset = 0;
        [y, x] = map2pix(refmat, easting-pixel_offset, northing+pixel_offset);
        [e,n]=pix2map(refmat,[floor(y) ceil(y)], [floor(x) ceil(x)]);
        val = imread(scene,'PixelRegion',{[floor(y) ceil(y)], [floor(x) ceil(x)]});
        val = double(val)./100;
        vals = max(0.05,interp2(e, n, val,easting,northing));
        dates = [];
        return
        
end

wildcard = strcat(fullfile(lsdir,folder,filesep), fbase,'*',fend,'.',ftype);

scenes = dir(wildcard);

names = extractfield(scenes, 'name')';

if iscellstr(names)
    names = char(names);
end

prs = str2num(names(:,5:9));
yrs = str2num(names(:, 10:13));
dys = str2num(names(:, 14:16));

imgs_idx = find(prs==pathrow);
scenes = scenes(imgs_idx);
dys = dys(imgs_idx);
yrs = yrs(imgs_idx);

obs = NaN(numel(scenes),3);
dates = NaT(numel(scenes),1);

%% Landsat pixels
for i=1:numel(scenes)
    dates(i) = datetime([yrs(i) 01 01]) + days(dys(i));
    fname = scenes(i).name;
    scene = fullfile(lsdir,folder,fname);
    info = geotiffinfo(scene);
    refmat = info.RefMatrix;
    
    [y, x] = map2pix(refmat, easting+15, northing-15);
    vals = imread(scene,'PixelRegion',{[floor(y) ceil(y)], [floor(x) ceil(x)]});
    vals = double(vals).*scale;
    [cfresult, conf] = is_cloudy(fname,y,x);
    cfresult = double(cfresult);
    conf = double(conf);
    [e,n]=pix2map(refmat,[floor(y) ceil(y)], [floor(x) ceil(x)]);
    e=e-15;
    n=n+15;
    if range(cfresult(:))==0 % all pixels have same cloud value
        px = interp2(e, n, vals,easting,northing);
        cfpx = mode(cfresult(:));
        conf = mean(conf(:));
    else % some pixel(s) with different cloud values
        if any(cfresult(:)==0) % at least 1 good pixel
            pts = [min(e) min(n); min(e) max(n); max(e) min(n);max(e) max(n)];
            pts = pts(cfresult==0,:);
            vals = vals(cfresult==0);
            
            pts(:,1) = pts(:,1) - easting;
            pts(:,2) = pts(:,2) - northing;
            dist = hypot(pts(:,1), pts(:,2));
            wt = 1./dist;
            px = sum(wt.*vals)./sum(wt);
            cfpx = 0.5;
            conf = mean(conf(:));
        else % no good pixels, but pixels with different values
            px = interp2(e, n, vals, easting, northing);
            cfpx = mode(cfresult(:));
            conf = mean(conf(:));
        end
    end
    obs(i,:) = [px cfpx conf];
end
land = obs(:,2)==0;
idw = obs(:,2)==0.5;
water = obs(:,2)==1;
shadow = obs(:,2)==2;
snow = obs(:,2)==3;
cloud = obs(:,2)==4;
forced = obs(:,2)==5;
good = land | idw | forced;

t = datenum(dates); % days
tq = min(t):max(t);

[xData, yData, wOut] = prepareCurveData(t(good), obs(good,1), 1-obs(good,2)./10);
dx = diff(xData);
dy = diff(yData);
spikes = abs([dy;0]./[dx;0])>spike_thresh & abs([0;dy]./[0;dx])>spike_thresh;
f = fit(xData(~spikes), yData(~spikes), interpmethod, 'Weights', wOut(~spikes));
xq = feval(f,tq);

ts = average_monthly_ts(tq, xq);
vals = ts{2};
dates = ts{1};
end

