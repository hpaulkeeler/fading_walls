function [sigMeanReal, sigMeanImag, sigVarReal, sigVarImag] = funSignalStats(sigInput)
% funSignalStats - Compute mean and variance of complex signal components
%
% Calculates basic statistics (mean and variance) separately for the real
% and imaginary parts of a complex electromagnetic signal.
%
% INPUT:
%   sigInput - Complex signal vector or matrix
%
% OUTPUTS:
%   sigMeanReal - Mean of real part: mean(Re{sigInput})
%   sigMeanImag - Mean of imaginary part: mean(Im{sigInput})
%   sigVarReal  - Variance of real part: var(Re{sigInput})
%   sigVarImag  - Variance of imaginary part: var(Im{sigInput})
%
% EXAMPLE:
%   S = randn(1000,1) + 1i*randn(1000,1);  % Complex Gaussian noise
%   [meanRe, meanIm, varRe, varIm] = funSignalStats(S);
%   fprintf('Real: mean=%.3f, var=%.3f\n', meanRe, varRe);
%   fprintf('Imag: mean=%.3f, var=%.3f\n', meanIm, varIm);
%
% NOTE:
%   For Rayleigh fading, Re{S} and Im{S} should both be zero-mean
%   Gaussian with equal variance

% Separate into real and imaginary components
sigReal = real(sigInput);
sigImag = imag(sigInput);

% Calculate means
sigMeanReal = mean(sigReal);
sigMeanImag = mean(sigImag);

% Calculate variances
sigVarReal = var(sigReal);
sigVarImag = var(sigImag);

end
