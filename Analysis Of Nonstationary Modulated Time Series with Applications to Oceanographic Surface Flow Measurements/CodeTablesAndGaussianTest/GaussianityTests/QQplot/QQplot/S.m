function [ out ] = S(params, N )
%This function returns the theoretical spectral density for the given
%parameters of the OU and Matérn, at Fourier frequencies.
%Parameters:
%params [1x6]:  parameters for the OU+Matérn. Note that the second
%               parameter, params(2), is the shift for the coriolis 
%               frequencies (these are already accounted for in the kernel).
%               Ask me about this point if you don't understand.
%Output:
%ESF [1xN]: the expected periodogram at Fourier frequencies (from -pi to
%pi).
out = S_((1-(0:N-1)/N), params);
end

