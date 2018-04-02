function [ output ] = kernel( gamma, Delta, T )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if Delta ~= 0
    output = 1/T*exp(1i*(gamma+Delta/(2*(T-1))) .*(0:T-1)) .* sin(Delta/(2*(T-1))*(0:T-1).*(T-(0:T-1)))./sin(Delta/(2*(T-1))*(0:T-1));
else
    output = (1-(0:T-1)/T).*exp(1i*gamma*(0:T-1));
end
output(1)=1;
end

