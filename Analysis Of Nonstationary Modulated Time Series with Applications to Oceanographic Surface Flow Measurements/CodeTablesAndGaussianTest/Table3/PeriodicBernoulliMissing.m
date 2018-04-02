function [ g ] = PeriodicBernoulliMissing(N, nbSamples, m, A, T)
%Function description:
%Returns
tempSeq = rand(N,nbSamples);
p = m+A/2*cos(2*pi/T*(1:N)');
g = zeros(N, nbSamples);
% for i=1:nbSamples
%     g(:,i) = (tempSeq(:,i) >= p);
% end
g = (tempSeq >= p);
end