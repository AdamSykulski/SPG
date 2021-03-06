%%In this script we test for the Gaussianity of the drifter velocities in
%%our data set of equatorial drifters.
clear all; close all;

%Add this file which has the Shapiro-Wilk test
addpath 'swtest'

%Load the drifters velocities
load blurreddrifters
load drifterulysses

%Parameters for the script
%'ks' (Kolmogorov-Smirnov) or 'sw' (Shapiro-Wilk)
test = 'ks';         
drifter_id = 1:200;
%Wether we differentiate the velocities
differentiate_velocities = 1;   
%Frequencies component above that one are set to 0 through a low-pass
%filter
highest_freq = 1;
%0 to test the longitudinal velocities, 1 to test
%the latitudinal velocities.
coordinate_test = 1;    

%To store the test results. The format is the following. Each row
%corresponds to one drifter velocity. Each row columns are as
%drifter_id, test at lag lag0 (see below), test at lags in distances.
results = zeros(length(drifter_id), 2);

%Do not change this
if strcmp(test, 'ks')
    test_function = @(x)kstest(x);
elseif strcmp(test, 'sw')
    test_function = @(x) swtest(x, 0.05);
end

%We compute the averaged autocovariance sequence
avg_ac = 0;     %initialization to any value;
for i=1:length(drifter_id)
    id = drifter_id(i);
    %We analyse drifters data one by one
    if id == 201
        X_ = drifterulysses.cv(1:852);
    else
        X_ = blurreddrifters.cv{id};
    end
    %Work on the filtered differentiated time series
    X_ = data_prepare(X_, differentiate_velocities, 1-highest_freq);
    
    %Sample autocovariance sequence
    ac = sample_autocorr(X_);
    %update the averaged autocorrelation sequence
    avg_ac = update_average_ac( avg_ac, ac, i);
end


%Find the first lag for which the autocovariance function takes a value
%smaller than the upper bound of the 95% confidence interval of the
%biased sample autocorrelation of a white noise process at same lag.
N = length(avg_ac);
lag0 = find(abs(avg_ac) > 1.96/sqrt(N), 1, 'last')+1;
lag0 = lag0/2;
disp(['Selected decorrelation lag :' num2str(lag0)]);
%Plot the average autocorrelation sequence
figure('name', 'Average autocorrelation sequence');
plot(abs(avg_ac));
hold on
line([0 N], [1.96/sqrt(N) 1.96/sqrt(N)], 'Color', 'red');

%Use that lag to make the tests
for i = 1:length(drifter_id)
    id = drifter_id(i);
    results(i, 1) = id;
    if id == 201
        X = drifterulysses.cv(1:852);
    else
        X = blurreddrifters.cv{id};
    end
    %Select one coordinate
    if coordinate_test
        X = -1i * X;
    end
    X = real(X);
    %Work on the low-pass filtered differentiated time series
    X = data_prepare(X, differentiate_velocities, 1-highest_freq);
    Y = X(1:lag0:end)/std(X(1:lag0:end));
    Y = Y(1:min(4000,end));
    h = test_function(Y);
    results(i, 2) = h;
end

%Compute the rejection rate
rejection_rate = sum(results(:,2:end))/length(drifter_id)*100;
rej_rate_lag0 = num2str(rejection_rate(1));
disp(['Rejection rate using selected subsampling lag: ' rej_rate_lag0]);