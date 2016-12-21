function [cmap] = rgbmap(varargin)
% RGBMAP creates color maps from any list of colors given as their common
% names.  Include a scalar anywhere in the list to specify the total number
% of levels in the color map.
% 
% This function requires the rgb function found on the Mathworks File
% Exchange site here: http://www.mathworks.com/matlabcentral/fileexchange/46872
% 
%
%% Syntax
% 
%  cmap = rgbmap('first color name','second color name')
%  cmap = rgbmap('first color name','second color name',...,'nth color name')
%  cmap = rgbmap(...,levels)
%  rgbmap(...)
%
%% Description 
% 
% |cmap = rgbmap('first color name','second color name')| creates an RGB
% color map |cmap| from some first color to a second color. 
% 
% |cmap = rgbmap('first color name','second color name',...,'nth color
% name')| creates a color map linearly scaled between any number of colors.
% 
% |cmap = rgbmap(...,M)| specifies the approximate number of levels |M| of the M x 3
% output colormap. Default value is 256. 
% 
% |rgbmap(...)| sets the color map without creating an array in the
% workspace. 
%
%
%% Example 1: Create a nice color map for some scattered data
% 
% x = 1:50; 
% y = cos(x*pi/50); 
% colors = rgbmap('blue','taupe','silver',50) ;
% scatter(x,y,30,colors,'filled')
%
%% Example 2: Get a 12-color map matrix from red to blue
% 
% cmap = rgbmap('red','blue',12)
% 
% h = surf(peaks);
% colorbar
% colormap(cmap)
% shading interp
% set(h,'edgecolor','k','edgealpha',.2)
% axis tight
% 
%% Example 3: Plot blue to white to red
% 
% h = surf(peaks);
% colorbar
% rgbmap('blue','white','red')
% shading interp
% set(h,'edgecolor','k','edgealpha',.2)
% caxis([-5 5])
% axis tight
% 
%% Author Info
% This function was written by Chad A. Greene on June 5, 2014. 
% 
% See also rgb and colormap. 
% 

%% Input check: 

assert(nargin>1,'Must have at least two inputs as strings.')
assert(exist('rgb','file')==2,'Cannot find the rgb function. Make sure Matlab knows where it is. If you don''t know where it is, you can download it here: http://www.mathworks.com/matlabcentral/fileexchange/46872')

%% Create color map: 

levels = 256; % default 
for k = 1:length(varargin)
    if isnumeric(varargin{k})
        levels = varargin{k};
        varargin(k)=[]; 
    end
end

for k = 1:length(varargin)
    try
        cp(k,:) = rgb(varargin{k}); 
    catch err
        error('MATLAB:rgbmap:rgb',['Cannot find an rgb value for the color ',varargin{k},'.'])
    end
end

numColors = length(varargin); 
levPerCol = floor(levels/(numColors - 1)); 

cmap=[]; 
for k = 1:numColors-1
    cmap(length(cmap)+1:length(cmap)+levPerCol,:) = [linspace(cp(k,1),cp(k+1,1),levPerCol)' linspace(cp(k,2),cp(k+1,2),levPerCol)' linspace(cp(k,3),cp(k+1,3),levPerCol)'];
end

% The next two while loops are a somewhat clunky fix in case the number of
% colors does not equal the number of levels specified by the user. If that
% happens, the last line(s) will be deleted or repeated until the correct
% number of colors is obtained. 
while length(cmap(:,1))>levels
    cmap(end,:)=[]; 
end
while length(cmap(:,1))<levels
    cmap(end+1,:)=cmap(end,:);
end

if nargout==0
    colormap(cmap)
    clear cmap
end

end

