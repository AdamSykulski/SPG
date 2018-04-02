function [acv] = complexouacvs(Q,N,Delta)
% computes the autocovariance for a complex OU of length N, spacing Delta,
% amplitude Q(1), frequency Q(2) and damping Q(3)
acv = Q(1)^2/abs(2*Q(3))*exp(1i*Q(2)*Delta*(0:N-1)).*exp(-abs(Q(3)*Delta*(0:N-1))); % acvs from lag 0 to N-1 (not 1 to N)