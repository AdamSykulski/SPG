function [ coriolisKernel ] = coriolis_kernel( corfreqs )
%This function returns the kernel that is used to compute the expected
%periodogram for the given sequence of coriolis frequencies.
%Parameters:
%corfreqs [1xN]: sequence of coriolis frequencies, with length N.
%Output:
%coriolisKernel [1xN]: the kernel that corresponds to the sequence of
%coriolis frequencies passed to the function.
N = length(corfreqs);
cumulative_rotation = exp(1i*cumsum(corfreqs));
coriolisKernel = zeros(1, N);
for k=1:N
    coriolisKernel(k) = 1/N*sum(cumulative_rotation(k:N).*conj(cumulative_rotation(1:N-k+1)));
end
end

