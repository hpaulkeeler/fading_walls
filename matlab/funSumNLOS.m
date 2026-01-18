function [nlosLoc, nlosPhase, nlosLocOpt, nlosPhaseOpt] = ...
    funSumNLOS(a, b, x, y, betaPath, waveNumb, kappaAbsorp, indexFirst, indexLast, numbOpt)
% funSumNLOS - Compute non-line-of-sight (NLOS) signal using method of images
%
% Calculates the sum of reflected electromagnetic rays between two parallel
% walls using the method of images. Returns signals for both random location
% and random phase models.
%
% INPUTS:
%   a, b          - Distances from origin to right and left walls
%   x, y          - Transmitter coordinates (can be vectors)
%   betaPath      - Path loss exponent β (signal decays as r^(-β/2))
%   waveNumb      - Wave number k = 2π/λ
%   kappaAbsorp   - Absorption coefficient κ ∈ [0,1] (power fraction absorbed)
%   indexFirst    - First reflection order to include (typically 1)
%   indexLast     - Last reflection order to include (e.g., 100)
%   numbOpt       - (Optional) Number of optimized phases (-1 for none)
%
% OUTPUTS:
%   nlosLoc       - NLOS signal for random location model
%   nlosPhase     - NLOS signal for random phase model
%   nlosLocOpt    - NLOS signal for optimized random location model
%   nlosPhaseOpt  - NLOS signal for optimized random phase model
%
% METHOD:
%   Uses method of images to compute infinite series of reflections.
%   The n-th reflection alternates between right and left walls:
%     - Right wall reflections: even indices
%     - Left wall reflections: odd indices
%   
%   Each reflection adds phase shift of π (inverts signal) and applies
%   absorption factor √κ. The series is truncated when terms become
%   negligible (< 1e-8).
%
% REFERENCE:
%   See manuscript "Signals reflecting off intelligent walls"
%   Section V: Two Parallel Walls

% Validate inputs
if ~isequal(size(x), size(y))
    error('x and y must have the same dimensions.');
end

% Reshape into column vectors
x = x(:);
y = y(:);

numbSim = length(x);                    % Number of locations to evaluate
numbRef = indexLast - indexFirst + 1;   % Number of reflections

%% ========================================================================
%% DISTANCE FUNCTIONS FOR METHOD OF IMAGES
%% ========================================================================

% Distance traveled by n-th reflection depends on whether ray last bounced
% off right or left wall, and whether total number of bounces is even/odd
%
% Pattern cycles through 4 cases:
%   1. Right-even: ray bounces even times, last bounce on right wall
%   2. Left-even:  ray bounces even times, last bounce on left wall  
%   3. Left-odd:   ray bounces odd times, last bounce on left wall
%   4. Right-odd:  ray bounces odd times, last bounce on right wall

% Right-even: 2i(a+b) + 2a - x
funRadRE = @(ii, x0, y0) sqrt(y0.^2 + (2*(ii-1)*(a+b) + 2*a - x0).^2);

% Left-even: 2i(a+b) + 2b + x
funRadLE = @(ii, x0, y0) sqrt(y0.^2 + (2*(ii-1)*(a+b) + 2*b + x0).^2);

% Right-odd: 2(k+1)(a+b) + x
funRadRO = @(kk, x0, y0) sqrt(y0.^2 + (2*kk*(a+b) + x0).^2);

% Left-odd: 2(k+1)(a+b) - x
funRadLO = @(kk, x0, y0) sqrt(y0.^2 + (2*kk*(a+b) - x0).^2);

%% ========================================================================
%% SUMMATION TERM FUNCTION
%% ========================================================================

% Calculate i-th term in the infinite sum for NLOS contribution
% Each reflection applies:
%   - Absorption factor: (-√κ)^i
%   - Path loss: r^(-β/2)
%   - Phase shift: exp(j*k*r)
funSumTerm = @(ii, kappa0, phase0, rad0) ...
    ((-sign(kappa0)*sqrt(abs(kappa0))).^ii) .* ...
    rad0.^(-betaPath/2) .* exp(1i*phase0);

%% ========================================================================
%% RANDOM PHASES FOR RANDOM PHASE MODEL
%% ========================================================================

% Generate uniformly distributed random phases on [0, 2π]
phaseRand = 2*pi * rand(numbSim, numbRef);

%% ========================================================================
%% INITIALIZE OUTPUT ARRAYS
%% ========================================================================

nlosLoc = zeros(numbSim, 1);        % Random location model
nlosPhase = zeros(numbSim, 1);      % Random phase model
nlosLocOpt = zeros(numbSim, 1);     % Optimized random location
nlosPhaseOpt = zeros(numbSim, 1);   % Optimized random phase

%% ========================================================================
%% MAIN SUMMATION LOOP
%% ========================================================================

% Determine which distance function to use (cycles through 4 cases)
numbRad = mod((0:indexLast-1), 4) + 1;

% Counters for indexing
ii_n = 1;  % Increments every 4 reflections
jj_n = 1;  % Increments every 2 reflections

for kk = indexFirst:indexLast
    
    % Select distance function based on reflection pattern
    switch numbRad(kk)
        case 1  % Right-even
            funRadAny = @(ii, x0, y0) funRadRE(ii, x0, y0);
        case 2  % Left-even
            funRadAny = @(ii, x0, y0) funRadLE(ii, x0, y0);
        case 3  % Left-odd
            funRadAny = @(ii, x0, y0) funRadLO(ii, x0, y0);
        case 4  % Right-odd
            funRadAny = @(ii, x0, y0) funRadRO(ii, x0, y0);
    end
    
    %% Random Location Model
    % Use actual ray distance with wave-based phase
    radLoc = funRadAny(ii_n, x, y);     % Distance traveled
    phaseLoc = waveNumb * radLoc;       % Phase = k*r
    termLoc = funSumTerm(jj_n, kappaAbsorp, phaseLoc, radLoc);
    nlosLoc = nlosLoc + termLoc;
    
    %% Random Phase Model
    % Use actual ray distance but with uniformly random phases
    radPhase = radLoc;                  % Same distance as location model
    phasePhase = phaseRand(:, kk);      % Uniformly distributed phases
    termPhase = funSumTerm(jj_n, kappaAbsorp, phasePhase, radPhase);
    nlosPhase = nlosPhase + termPhase;
    
    %% Optimized Models (if requested)
    if exist('numbOpt', 'var') && numbOpt > 1
        % Placeholder for optimized phase calculations
        % (Implementation depends on optimization strategy)
        
        if kk <= numbOpt
            % Update with optimized terms
            nlosLocOpt = nlosLocOpt + termLoc;
            nlosPhaseOpt = nlosPhaseOpt + termPhase;
        else
            % Update with non-optimized terms
            nlosLocOpt = nlosLocOpt + termLoc;
            nlosPhaseOpt = nlosPhaseOpt + termPhase;
        end
    end
    
    %% Early Termination Check
    % Stop if contribution becomes negligible
    termAbsMax = max(abs(termLoc));
    tolTerm = 1e-8;
    
    if termAbsMax < tolTerm
        break;  % Remaining terms will be even smaller
    end
    
    %% Update Counters
    ii_n = ii_n + (mod(kk, 4) == 0);  % Increment every 4 reflections
    jj_n = jj_n + (mod(kk, 2) == 0);  % Increment every 2 reflections
    
end

end
