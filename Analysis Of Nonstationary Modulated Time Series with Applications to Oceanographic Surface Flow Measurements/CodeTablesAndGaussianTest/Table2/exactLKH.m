function [ output ] = exactLKH( a, sigma, gamma, Delta, data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
T = length(data);
output1 = T*log(1/(pi*sigma^2));
output2 = 1/sigma^2*(abs(data(1))^2+sum(abs(data(2:end)-a*exp(1i*linspace(gamma-Delta/2, gamma+Delta/2, T-1)').*data(1:T-1)).^2));
output = -(output1-output2);
end

