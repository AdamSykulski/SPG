function [ ker ] = kernel(g)
%Function description:
%Computes the kernel corresponding to the modulating sequence that is
%passed as a parameter.
%
%Parameters:
%g [Nx1] g modulating sequence.
N = length(g);
% ker = zeros(1,N);
% for k=1:N
%     ker(k) = 1/N*sum(g(1:N-k+1).*g(k:N));
% end
f = 1/(2*N)*abs(fft(g,2*N)).^2;
ker2 = ifft(2*f);
ker = ker2(1:N)';
end

