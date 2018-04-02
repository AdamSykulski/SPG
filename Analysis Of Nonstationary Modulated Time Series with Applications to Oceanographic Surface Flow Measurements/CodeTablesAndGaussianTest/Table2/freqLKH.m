function [ output ] = freqLKH( a, sigma, gamma, Delta, periodogram )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
T = length(periodogram);
ker = kernel(gamma, Delta, T);
c = AR1C.autocov_seq(a, sigma, T-1, false);
c_ = ker .* c;
S_ = 2*real(fft(c_))-c_(1);
output = sum(log(S_)+periodogram'./S_);
end

