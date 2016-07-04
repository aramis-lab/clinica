function map = spectral( m );

%Black-purple-blue-green-yellow-red-white color map.
%
% Usage: map = spectral( m );
%
% m = number of colours; default is the length of the current colormap.
%
% map = m x 3 matrix containing a "spectral" colormap.
% 
% E.g. colormap(spectral) resets the colormap of the current figure.         

if nargin < 1, m = size(get(gcf,'colormap'),1); end
base = [
  0.2000 0.2000 0.2000
  0.4667 0.0000 0.5333
  0.5333 0.0000 0.6000
  0.0000 0.0000 0.6667
  0.0000 0.0000 0.8667
  0.0000 0.4667 0.8667
  0.0000 0.6000 0.8667
  0.0000 0.6667 0.6667
  0.0000 0.6667 0.5333
  0.0000 0.6000 0.0000
  0.0000 0.7333 0.0000
  0.0000 0.8667 0.0000
  0.0000 1.0000 0.0000
  0.7333 1.0000 0.0000
  0.9333 0.9333 0.0000
  1.0000 0.8000 0.0000
  1.0000 0.6000 0.0000
  1.0000 0.0000 0.0000
  0.8667 0.0000 0.0000
  0.8000 0.0000 0.0000
  0.8000 0.8000 0.8000
];
n = length(base);
X0 = linspace (1, n, m);
map = interp1(1:n,base,X0);

return
