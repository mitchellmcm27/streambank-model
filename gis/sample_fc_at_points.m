function [vals] = sample_fc_at_points(points, pathrow)


prstring = num2str(pathrow);
    scene = strcat(fullfile('data','forest',filesep),...
        'p0',prstring(1:2),'r',prstring(3:5),'_TC_2005.tif');
    info = geotiffinfo(scene);
    refmat = info.RefMatrix;
    %pixel_offset = -15;
    pixel_offset = 0;
    
vals = nan(size(points,1),1);
for i=1:numel(vals)
    
    point = points(i,:);
    easting = point(1);
    northing = point(2);
    
    [y, x] = map2pix(refmat, easting-pixel_offset, northing+pixel_offset);
    [e,n]=pix2map(refmat,[floor(y) ceil(y)], [floor(x) ceil(x)]);
    val = imread(scene,'PixelRegion',{[floor(y) ceil(y)], [floor(x) ceil(x)]});
    val = double(val)./100;
    vals(i) = max(0.05,interp2(e, n, val,easting,northing));
end

vals

end

