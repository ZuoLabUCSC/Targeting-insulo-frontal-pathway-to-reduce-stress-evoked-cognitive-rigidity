function cmap = green_white_magenta(n)
%GREEN_WHITE_MAGENTA Creates a green-to-white-to-magenta colormap
% usage example:
% colormap(green_white_magenta(256));
% colorbar;
%
%   cmap = GREEN_WHITE_MAGENTA(n) returns an n-by-3 colormap array.
%   If n is not specified, it defaults to 256.

    if nargin < 1
        n = 256; % Default number of colors
    end

    % Define the key colors
    green = [0, 1, 0];
    white = [1, 1, 1];
    magenta = [1, 0, 1];

    % Interpolate between green and white
    n1 = round(n/2);
    cmap1 = [linspace(green(1), white(1), n1)', ...
             linspace(green(2), white(2), n1)', ...
             linspace(green(3), white(3), n1)'];

    % Interpolate between white and magenta
    n2 = n - n1;
    cmap2 = [linspace(white(1), magenta(1), n2)', ...
             linspace(white(2), magenta(2), n2)', ...
             linspace(white(3), magenta(3), n2)'];

    % Combine the two parts
    cmap = [cmap1; cmap2];
end