%  SHARED UTILITY: red-blue diverging colormap
function cmap = redblue_cmap(n)
if nargin < 1, n = 256; end
r = [linspace(0,1,n/2), ones(1,n/2)];
b = [ones(1,n/2), linspace(1,0,n/2)];
g = [linspace(0,1,n/2), linspace(1,0,n/2)];
cmap = [r' g' b'];
end
