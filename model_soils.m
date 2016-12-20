function [val] = model_soils(point)
soilsdir = fullfile('K:','GIS','model','input', 'soils');

easting = point(1);
northing = point(2);

bdfile = strcat(fullfile(soilsdir,'AVG_BD_UTM16N.tif'));
info = geotiffinfo(bdfile);
refmat = info.RefMatrix;
[y, x] = map2pix(refmat, easting, northing);
vals = imread(bdfile, 'PixelRegion',{[floor(y) ceil(y)], [floor(x) ceil(x)]});
vals = double(vals);
vals(vals<0) = NaN;
if ~any(isnan(vals))
    [e,n]=pix2map(refmat,[floor(y) ceil(y)], [floor(x) ceil(x)]);
    val = interp2(e, n, vals, easting, northing);
else
    step=0;
    while all(all(isnan(vals)))
        step=step+1;
        vals = imread(bdfile, 'PixelRegion',{[floor(y)-step ceil(y)+step], [floor(x)-step ceil(x)+step]});
        vals = double(vals);
        vals(vals<0) = NaN;
    end
    val = mode(mode(vals));
    
    
end
end