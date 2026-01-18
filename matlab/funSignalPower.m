function [sigPow, sigReal, sigImag] = funSignalPower(sigInput, expPower)
% funSignalPower - Extract power, real, and imaginary parts of signal
%
% Computes the power (magnitude raised to an exponent) and separates the
% real and imaginary components of a complex electromagnetic signal.
%
% INPUTS:
%   sigInput  - Complex signal vector or matrix
%   expPower  - (Optional) Exponent for power calculation (default: 2)
%               Power = |signal|^expPower
%
% OUTPUTS:
%   sigPow   - Signal power: |sigInput|^expPower
%   sigReal  - Real part of signal: Re{sigInput}
%   sigImag  - Imaginary part of signal: Im{sigInput}
%
% TYPICAL USE:
%   expPower = 2 gives P = |S|^2 (signal power)
%   expPower = 1 gives magnitude |S| (signal amplitude)
%
% EXAMPLE:
%   S = [1+2i, 3-i, 2+0i];
%   [P, Re_S, Im_S] = funSignalPower(S, 2);
%   % Result: P = [5, 10, 4], Re_S = [1, 3, 2], Im_S = [2, -1, 0]

% Set default power exponent
if nargin == 1
    expPower = 2;
end

% Extract real and imaginary parts
sigReal = real(sigInput);
sigImag = imag(sigInput);

% Calculate power: P = |S|^expPower
sigPow = abs(sigInput).^expPower;

end
