function [new_average_ac,  new_n ] = update_average_ac( average_ac, new_ac, n )
%This function updates an average of the autocorrelation sequences of
%independent realizations of time series. The updated average length is the
%minimum of the previous average length and of the new autocorrelation
%sequence.
if n==1
    N = length(new_ac);
    average_ac = zeros(N,1);
else
    N = min(length(average_ac), length(new_ac));
    average_ac = average_ac(1:N);
end
new_ac = new_ac(1:N);
new_average_ac = 1/n*new_ac + (n-1)/n*average_ac;
new_n = n+1;
end

