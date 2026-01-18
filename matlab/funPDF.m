function [f, xi] = funPDF(X, booleMethod)
% funPDF - Estimate probability density function from data
%
% Estimates the PDF of data using either histogram binning with smoothing
% or kernel density estimation.
%
% INPUTS:
%   X           - Data vector from which to estimate PDF
%   booleMethod - (Optional) PDF estimation method:
%                 1: Histogram with 100 bins + 3-point moving average smoothing
%                 0 or omitted: Kernel density estimation (requires Statistics Toolbox)
%
% OUTPUTS:
%   f  - Estimated probability density values
%   xi - Corresponding x-values where density is evaluated
%
% METHOD 1 (booleMethod = 1):
%   - Creates 100-bin histogram
%   - Normalizes to PDF (area = 1)
%   - Applies 3-point moving average for smoothing
%   - Faster but may miss fine details
%
% METHOD 0 (booleMethod = 0 or omitted):
%   - Uses ksdensity for kernel density estimation
%   - Smoother result but slower
%   - Requires Statistics and Machine Learning Toolbox
%
% EXAMPLE:
%   X = randn(10000, 1);  % Normal random data
%   [f, xi] = funPDF(X, 1);  % Histogram method
%   plot(xi, f);
%
% NOTE:
%   Histogram method is recommended for large datasets (N > 10000)
%   Kernel method is better for smaller datasets or when smoothness matters

if (nargin == 2) && (booleMethod == 1)
    % METHOD 1: Histogram with smoothing
    
    % Create histogram with 100 bins, normalized as PDF
    [pdfNEmp, vecX] = histcounts(X, 100, 'Normalization', 'pdf');
    
    % Apply 3-point moving average to smooth histogram
    % This reduces noise while preserving overall shape
    pdfNEmp = movmean(pdfNEmp, 3);
    
    % Output
    f = pdfNEmp;
    xi = vecX(1:end-1);  % Bin centers (left edges)
    
else
    % METHOD 0: Kernel density estimation
    % Note: Requires Statistics and Machine Learning Toolbox
    [f, xi] = ksdensity(X);
end

end
