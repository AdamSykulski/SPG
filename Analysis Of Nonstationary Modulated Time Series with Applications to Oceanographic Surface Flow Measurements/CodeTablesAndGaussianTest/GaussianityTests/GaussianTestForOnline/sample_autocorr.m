function [ acr ] = sample_autocorr( X )
%This function returns the biased autocovariance sequence, based on a FFT,
%making it a N*log(N) computational order.
N = length(X);
P = 1/(2*N)*abs(fft(X, 2*N)).^2;
acr = ifft(2*P);
acr = acr(1:N);
acr = acr / acr(1);
end