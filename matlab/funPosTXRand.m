function [xRand, yRand] = funPosTXRand(xTX, yTX, xDeltaTX, yDeltaTX, choiceVaryXY, numbRand)
% funPosTXRand - Generate random transmitter locations for empirical sampling
%
% Randomly samples transmitter locations uniformly within specified regions.
% Used to generate empirical probability density functions (PDFs) from
% Monte Carlo simulation.
%
% INPUTS:
%   xTX          - Vector of x-coordinates for transmitter centers [x1, x2, x3]
%   yTX          - Vector of y-coordinates for transmitter centers [y1, y2, y3]
%   xDeltaTX     - Vector of x sampling widths [Δx1, Δx2, Δx3]
%   yDeltaTX     - Vector of y sampling widths [Δy1, Δy2, Δy3]
%   choiceVaryXY - Sampling mode:
%                  1: vary x only (use xTX(1), xDeltaTX(1))
%                  2: vary y only (use yTX(2), yDeltaTX(2))
%                  3: vary both x and y (use xTX(3), yTX(3), etc.)
%   numbRand     - Number of random samples to generate
%
% OUTPUTS:
%   xRand - Column vector of random x-coordinates (length numbRand)
%   yRand - Column vector of random y-coordinates (length numbRand)
%
% SAMPLING REGIONS:
%   Mode 1 (vary x): xRand ~ Uniform[x - Δx/2, x + Δx/2], yRand = y (fixed)
%   Mode 2 (vary y): xRand = x (fixed), yRand ~ Uniform[y - Δy/2, y + Δy/2]
%   Mode 3 (vary x and y): Both coordinates uniformly random in rectangle
%
% EXAMPLE:
%   % Generate 10000 random locations varying x around x=0.25 ± 0.1
%   xTX = [0.25, 0, 0.25];
%   yTX = [0, 0.25, 0.25];
%   xDeltaTX = [0.2, 0, 0.2];
%   yDeltaTX = [0, 0.2, 0.2];
%   [xRand, yRand] = funPosTXRand(xTX, yTX, xDeltaTX, yDeltaTX, 1, 10000);
%   % Result: xRand ~ Uniform[0.15, 0.35], yRand = 0
%
% SEE ALSO:
%   funPosTXPlot - Generate deterministic grid of locations for plotting

% Extract parameters for selected mode
x = xTX(choiceVaryXY);
y = yTX(choiceVaryXY);
xDelta = xDeltaTX(choiceVaryXY);
yDelta = yDeltaTX(choiceVaryXY);

switch choiceVaryXY
    case 1
        % Vary x only, keep y fixed
        xMin = x - xDelta/2;
        xRand = xDelta * rand(numbRand, 1) + xMin;
        yRand = y * ones(numbRand, 1);
        
    case 2
        % Keep x fixed, vary y only
        yMin = y - yDelta/2;
        xRand = x * ones(numbRand, 1);
        yRand = yDelta * rand(numbRand, 1) + yMin;
        
    case 3
        % Vary both x and y (2D rectangular region)
        xMin = x - xDelta/2;
        yMin = y - yDelta/2;
        xRand = xDelta * rand(numbRand, 1) + xMin;
        yRand = yDelta * rand(numbRand, 1) + yMin;
        
    otherwise
        error('choiceVaryXY must be 1, 2, or 3');
end

end
