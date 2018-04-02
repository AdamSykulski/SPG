function [ res ] = S_( c, ker )
%Computes the expected periodogram for a finite sample of a modulated time
%series. The finite autocovariance sequence of the latent process is passed
%as a parameter, as well as the kernel corresponding to the modulating
%sequence.
%
%Parameters:
%c [Nx1] autocovariance sequence
%ker [Nx1] kernel sequence corresponding to the modulation sequence
c_ = c .* ker;
res = 2*real(fft(c_))-c_(1);
end

