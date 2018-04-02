function [stop] = add_frame( periodogram, theta, ker, frames, optim_values)
%UNTITLED37 Summary of this function goes here
%   Detailed explanation goes here
r = theta(1); 
sigma = theta(2);
N = length(periodogram);
expected = S_(acvsAR1(r, 0, sigma, N), ker);
hold off
plot(10*log10(periodogram));
hold on
plot(10*log10(expected));
drawnow 
i_frame = optim_values.iteration;
frames(i_frame+1) = getframe(gcf);
stop=false;
end

