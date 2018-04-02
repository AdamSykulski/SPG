function [ acvs ] = acvsAR1(r, theta, sigma, N )
%Function description:
%Returns the finite autocovariance sequence of a complex-valued AR(1)
%process with damping, frequency and noise amplitude parameters passed as
%arguments. Note that if the frequency parameter is zero then the function returns
%the autocovariance sequence of a real-valued AR(1) process.
%
%Parameters:
%r              float:      damping parameter
%theta          float:      frequency parameter
%sigma          float:      noise amplitude
%N              int:        time series length. The function returns the 
%                           autocovariance sequence for lags 0, ..., N-1.
acvs = sigma^2/(1-r^2)*r.^((0:N-1)).*exp(1i*theta*(0:N-1));
end