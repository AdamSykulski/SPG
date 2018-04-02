function [ params, fval, exit_flag, CPUtime, ker ] = fitOUandMatern( SZ, corfreqs, N, LB, UB, MF, ZEROF, start, minBound, maxBound, options)
%This function fits the OU and Matern model using the nonstationary
%likelihood.
%Parameters:
%SZ [1xN]: the periodogram
%corFreqs [1xN]: the sequence of coriolis frequencies
%N: the length of the sample
%LB [1x1]: lower bound for the frequency estimation
%UP [1x1]: upper bound for the frequency estimation
%MF
%ZEROF [1x1]: whether zero frequency should be included
%start [1x6]: starting parameters for the estimation
%minBound [1x6]: minimum bound for the parameters
%maxBound [1x6]: maximum bound for the parameters


%We start the time counter
tic;

%We start by obtaining the kernel which corresponds to the sequence of
%coriolis frequencies passed to the function.
ker = coriolis_kernel(corfreqs);

%We minimize the -likelihood
[xr, fval, exit_flag] = fminsearchbnd(@(x) lkh_(x, start, SZ, ker, N,LB,UB,MF,ZEROF), ones(1,6), minBound, maxBound,options);

%We obtain the estimated params by multiplying the starting parameters by
%the estimated factor.
params = xr .* start;

%We also return the CPU time
CPUtime = toc;
end

