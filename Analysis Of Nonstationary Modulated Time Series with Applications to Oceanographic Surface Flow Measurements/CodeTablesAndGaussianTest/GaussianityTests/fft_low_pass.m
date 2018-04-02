function [ Y ] = fft_low_pass( X, p )
%This function acts as a low pass filter. The signal X is Fourier
%transformed, then the 100*p percent higher frequencies are set to zero.
%The filtered signal is obtained by inverse Fourier transforming.
if p==1
    Y = X;
else
    N = length(X);
    freq_domain = fft(X);
    freq_domain = fftshift(freq_domain);
    freq_domain(1:floor(N/2*p)) = 0;
    freq_domain(ceil(N/2*(2-p)):end) = 0;
    Y = real(ifft(ifftshift(freq_domain)));
end
end

