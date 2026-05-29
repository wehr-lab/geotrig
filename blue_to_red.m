function c = blue_to_red(n)

if nargin < 1
    n = 256; %default
end

% Define the colors (using 0-1 RGB triplets)
blue = [0, 0, 1];
grey = [0.75, 0.75, 0.75];
red = [1, 0, 0];

% Create the blue-to-grey gradient (n/2 steps)
c1 = [linspace(blue(1), grey(1), n/2)', ...
      linspace(blue(2), grey(2), n/2)', ...
      linspace(blue(3), grey(3), n/2)'];

% Create the grey-to-red gradient (n/2 steps)
c2 = [linspace(grey(1), red(1), n/2)', ...
      linspace(grey(2), red(2), n/2)', ...
      linspace(grey(3), red(3), n/2)'];

% Combine them into one colormap
c = [c1; c2];

% Apply it
%colormap(my_cmap);
%colorbar;