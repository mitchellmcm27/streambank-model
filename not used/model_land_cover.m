function [lc] = model_land_cover(lat, lon, varargin)
old_dir = cd;
cd('K:\GIS\r\USGS NLCD\nlcd_2011_landcover_2011_edition_2014_10_10\');

if nargin == 3
    utmzone = varargin{1};
    [lat, lon] = utm2deg(lat, lon, utmzone);
end

info = geotiffinfo('nlcd2');
[x, y] = projfwd(info,lat,lon);
[row, col] = map2pix(info.RefMatrix,x,y);
row = row + 0.5;
col = col + 0.5;
lc = imread('nlcd2.tif',...
    'PixelRegion', {[round(row) round(row)], [round(col) round(col)]});

lc = double(lc);

i=1;
test = lc;
while test<31
    region = imread('nlcd2.tif',...
        'PixelRegion', {[floor(row-i) ceil(row+i)], [floor(col-i) ceil(col+i)]});
    test = mode(region(region>11));
    if isnan(test)
        test = lc;
    end
    i=i+1;
end
lc = test;
cd(old_dir);
end
