function [ res ] = lkh_(P, r, sigma, N, ker)
%pseudo-likelihood function.
%Returns the evaluation of the pseudo-likelihood for the periodogram
%and a fitted AR(1). The kernel corresponding to the modulating sequence
%is passed as a parameter.
%
%Parameters:
%P [Nx1]:       periodogram of the time series
%r int:         damping parameter of the fitted AR(1)
%sigma int:     variance of the noise of the fitted AR(1)
%N int:         size of the time series
%ker [Nx1]:     kernel
%
c = acvsAR1(r, 0, sigma, N );
S = S_(c, ker);
res = sum(log(S)+P./S);
end

