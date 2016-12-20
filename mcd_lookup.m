function [ mcd ] = mcd_lookup( nlcd )
% Convert USGS NLCD land cover types to MODIS MCD12Q1 biomes using the
% table 2 of Ganguly et al 2012, RSE

if nlcd==11 || nlcd==12 % open water or ice/snow
    mcd = 254;
elseif nlcd==31 || nlcd==32 % barren/rocky
    mcd = 255;
elseif nlcd==42 || nlcd==43 % evergreen or mixed
    mcd = 7; % evergreen needleleaf forest
    % could also be 5, evergreen broadleaf forest
elseif nlcd==41
    mcd = 6; % deciduous broadleaf forest
    % could also be 8, deciduous needleleaf forest
elseif nlcd==90 || nlcd==95 % woody wetlands
    mcd = 4; % savanna??? what
elseif nlcd==21 || nlcd==22 || nlcd==23 % developed
    mcd = 4; % savanna
elseif nlcd == 82 % cultivated crops
    mcd = 3; % broadleaf crops
elseif nlcd == 51 || nlcd == 52 % shrub, scrub
    mcd = 2; % shrubs
elseif nlcd == 71 || nlcd==72 || nlcd==73 || nlcd==74 || nlcd==81
    mcd = 1; % grasses/cereal crops
end


end

