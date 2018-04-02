function [ Y ] = data_prepare( X, differentiate, freq_cutoff )
%This function prepares the data. If differentiate is 1, then the time
%series is differentiated. freq_cutoff is between 0 and 1 and represents
%the percentage of frequencies set to zero.
Y = X;
if differentiate
    Y = diff(X);
end
Y = fft_low_pass(Y, freq_cutoff);
end

