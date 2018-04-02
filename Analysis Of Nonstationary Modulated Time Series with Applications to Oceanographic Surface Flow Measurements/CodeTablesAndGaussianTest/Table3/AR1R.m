function [ data ] = AR1R( r, sigma, N, nbSamples )
%Function description:
%Generates independent samples according to a real-valued AR1 model with
%given damping and noise amplitude parameters.
%
%Parameters:
%r          float:              damping parameter
%sigma      float:              noise variance
%N          int:                time series length
%nbSamples  int:                number of independent samples
%
%Output:
%data       float[nbSamplesxN]  Generated time series. Each line
%                               corresponds to a sample.
data = zeros(nbSamples, N);
%Initialization
data(:,1) = randn(nbSamples, 1)*sigma/sqrt((1-r^2));
noise = randn(nbSamples, N)*sigma;

%Generation of the data
for i=1:N-1
    data(:, i+1) = r*data(:,i) + noise(:, i+1);
end
end

