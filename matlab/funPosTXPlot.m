function [xPlot, yPlot] = funPosTXPlot(xTX, yTX, xDeltaTX, yDeltaTX, choiceVaryXY, numbPlot)
% funPosTXPlot - Generate deterministic grid of transmitter locations for plotting
%
% Creates a grid of evenly-spaced transmitter locations for plotting 
% deterministic power curves and surfaces.
%
% INPUTS:
%   xTX          - Vector of x-coordinates for transmitter centers [x1, x2, x3]
%   yTX          - Vector of y-coordinates for transmitter centers [y1, y2, y3]
%   xDeltaTX     - Vector of x plotting widths [Δx1, Δx2, Δx3]
%   yDeltaTX     - Vector of y plotting widths [Δy1, Δy2, Δy3]
%   choiceVaryXY - Plotting mode:
%                  1: vary x only (1D line plot)
%                  2: vary y only (1D line plot)
%                  3: vary both x and y (2D surface plot)
%   numbPlot     - Number of points to generate
%
% OUTPUTS:
%   xPlot - X-coordinates for plotting (column vector or matrix)
%   yPlot - Y-coordinates for plotting (column vector or matrix)
%
% BEHAVIOR BY MODE:
%   Mode 1 (vary x): Creates numbPlot points along x-axis, y fixed
%                    Output: xPlot, yPlot are column vectors
%
%   Mode 2 (vary y): Creates numbPlot points along y-axis, x fixed
%                    Output: xPlot, yPlot are column vectors
%
%   Mode 3 (vary x and y): Creates meshgrid with sqrt(numbPlot)*10 points
%                          in each dimension
%                          Output: xPlot, yPlot are matrices from meshgrid
%
% EXAMPLE:
%   % Generate 500 points along x for plotting P(x, y=0)
%   xTX = [0.25, 0, 0.25];
%   yTX = [0, 0.25, 0.25];
%   xDeltaTX = [0.2, 0, 0.2];
%   yDeltaTX = [0, 0.2, 0.2];
%   [xPlot, yPlot] = funPosTXPlot(xTX, yTX, xDeltaTX, yDeltaTX, 1, 500);
%   % Result: xPlot = linspace(0.15, 0.35, 500)', yPlot = zeros(500,1)
%
% SEE ALSO:
%   funPosTXRand - Generate random locations for empirical PDF estimation

% Extract parameters for selected mode
x = xTX(choiceVaryXY);
y = yTX(choiceVaryXY);
xDelta = xDeltaTX(choiceVaryXY);
yDelta = yDeltaTX(choiceVaryXY);

switch choiceVaryXY
    case 1
        % Vary x only, keep y fixed (1D plot)
        xMin = x - xDelta/2;
        xMax = x + xDelta/2;
        xPlot = linspace(xMin, xMax, numbPlot)';
        yPlot = y * ones(numbPlot, 1);
        
    case 2
        % Keep x fixed, vary y only (1D plot)
        yMin = y - yDelta/2;
        yMax = y + yDelta/2;
        xPlot = x * ones(numbPlot, 1);
        yPlot = linspace(yMin, yMax, numbPlot)';
        
    case 3
        % Vary both x and y (2D surface plot)
        xMin = x - xDelta/2;
        xMax = x + xDelta/2;
        yMin = y - yDelta/2;
        yMax = y + yDelta/2;
        
        % Create meshgrid with resolution proportional to numbPlot
        numbMeshX = round(10 * sqrt(numbPlot));
        numbMeshY = round(10 * sqrt(numbPlot));
        xMesh = linspace(xMin, xMax, numbMeshX);
        yMesh = linspace(yMin, yMax, numbMeshY);
        [xPlot, yPlot] = meshgrid(xMesh, yMesh);
        
    otherwise
        error('choiceVaryXY must be 1, 2, or 3');
end

end
