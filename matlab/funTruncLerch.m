function y = funTruncLerch(z, s, v)
% funTruncLerch - Compute the Lerch transcendent function (modified version)
%
% Calculates a truncated approximation of the Lerch transcendent function,
% which is defined as the infinite power series:
%
%   Φ(z, s, v) = Σ_{n=0}^∞ z^n / (v + n)^s
%
% This implementation starts from n=1 instead of n=0, effectively computing:
%
%   Φ_modified(z, s, v) = Σ_{n=1}^∞ z^n / (v + n)^s
%
% This is appropriate for computing the NLOS (non-line-of-sight) signal,
% where the n=0 term represents the direct LOS signal and is handled separately.
%
% INPUTS:
%   z - Complex parameter (typically involves absorption and wave number)
%   s - Power parameter (typically β/2, where β is path loss exponent)
%   v - Offset parameter (typically ±r/d, where r is distance and d = a+b)
%       Can be a scalar or vector
%
% OUTPUT:
%   y - Truncated Lerch function value(s), same size as v
%
% PARAMETERS:
%   numbFirst - First term index (set to 1 for NLOS calculation)
%   numbLast  - Last term index (400 provides good convergence)
%
% NOTES:
%   - Standard Lerch function starts at n=0; this version starts at n=1
%   - The n=0 term (1/v^s) must be subtracted separately if using standard
%     definition (see manuscript equation 29)
%   - For vectorized input v, output is reshaped to match input dimensions
%
% REFERENCE:
%   NIST Digital Library of Mathematical Functions
%   https://dlmf.nist.gov/25.14
%
% SEE ALSO:
%   For standard Lerch function (starting from n=0), use numbFirst=0

% Summation range
numbFirst = 1;      % Start from n=1 (excludes LOS term)
numbLast = 400;     % Sufficient terms for convergence

% Index vector
nn = numbFirst:numbLast;

if length(v) > 1
    % Vectorized calculation for multiple v values
    vInput = v;
    v = v(:);  % Ensure column vector
    
    % Create matrices for vectorized computation
    nnMatrix = repmat(nn, length(v), 1);        % Each row is [1,2,...,400]
    vMatrix = repmat(v, 1, length(nn));         % Each col is v values
    
    % Compute sum: Σ z^n / (v+n)^s
    y = sum((z.^nnMatrix) ./ (vMatrix + nnMatrix).^s, 2);
    
    % Reshape output to match input shape
    y = reshape(y, size(vInput));
    
else
    % Scalar calculation
    y = sum((z.^nn) ./ (v + nn).^s);
end

end
