%%In this script we test for the Gaussianity of the drifter velocities in
%%our data set of equatorial drifters.
clear all; close all;

%Add this file which has the Shapiro-Wilk test
addpath 'swtest'

%Load the drifters velocities
load ../blurreddrifters
load ../drifterulysses

%fig1 = figure('Name', 'X coordinate');

%Parameters for the script
%'ks' (Kolmogorov-Smirnov) or 'sw' (Shapiro-Wilk)
test = 'ks';         


drifter_id = 1:201;
distances = 1:200;
%ids for which we want a plot of the empirical cumulative distribution. Do
%not chose too many of them.
plot_ids = 1:5;
plot_distances = [20 50 100];

%To store the test results. The format is the following. Each row
%corresponds to one drifter velocity. Each row columns are as
%drifter_id, test at lag lag0 (see below), test at lags in distances.
results = zeros(length(drifter_id), length(distances)+2);

%Do not change this
if strcmp(test, 'ks')
    test_function = @(x)kstest(x);
elseif strcmp(test, 'sw')
    test_function = @(x) swtest(x, 0.05);
else
    %I get a warning message when using this one so I've left it aside for
    %now.
    test_function = @(x) adtest(x);
end

%Manage plot ids
if ~isempty(plot_ids)
    next_plot_id = plot_ids(1);
    n_plot_id = 1;
else
    next_plot_id = -1;
end

for i=1:length(drifter_id);
    %We analyse drifters data one by one
    id = drifter_id(i)
    results(i,1) = id;
    if id == 201
        X_ = drifterulysses.cv(1:852);
        X = real(drifterulysses.cv(1:852));
    else
        X_ = blurreddrifters.cv{id};
        X = real(blurreddrifters.cv{id});
    end
    %Work on the differentiated time series
    X = diff(X);
    X_ = diff(X_);
    X = fft_low_pass(X, 1-0.8/6);
    X_ = fft_low_pass(X_, 1-0.8/6);
    N = length(X);
    
    %Sample autocovariance sequence
    ac = sample_autocorr(X_);
    %Find the first lag for which the autocovariance function takes a value
    %smaller than the upper bound of the 95% confidence interval of the
    %biased sample autocorrelation of a white noise process at same lag.
    lag0 = find(abs(ac) < 1.96/sqrt(N), 1, 'first')
    distances_i = [lag0 distances];
    plot_distances_i = [lag0 plot_distances];
    next_plot_dist = lag0;
    %Plot the sample autocorrelation for the first one
    if i==1
        figure('name', 'Sample autocorrelation sequence');
        plot(abs(ac));
        hold on
        line([0 N], [1.96/sqrt(N) 1.96/sqrt(N)], 'Color', 'red');
    end
    
    if i==next_plot_id
        figure('name', ['Drifter #' num2str(id)]);
    end
    
    
    %Test for the different distances. distances_i(1) is lag0.
    subplot_nb = 1;
    for j=1:length(distances_i)
        d = distances_i(j);
        Y=X(1:d:end)/std(X(1:d:end));
        Y = Y(1:min(4000,end));
        h = test_function(Y);
        results(i, j+1) = h;
        if(i==next_plot_id && j==next_plot_dist)
            %Plot of the empirical cdf. Only for the 4 first lags.
            subplot(220+subplot_nb);
            fig_title= ['Lag: ' num2str(d) ' - subsample size: ' num2str(length(Y)) ' - Test rejected: ', num2str(h)];
            plot_empirical_cdf(Y, fig_title);
            subplot_nb = subplot_nb + 1;
            if subplot_nb>4
                next_plot_dist = -1;
            else
                next_plot_dist = plot_distances_i(subplot_nb);
            end
        end
    end
    
    %Update the plot id
    if i == next_plot_id
        if length(plot_ids) > n_plot_id
            n_plot_id = n_plot_id + 1;
            next_plot_id = plot_ids(n_plot_id);
        end
    end
end
rejection_rate = sum(results(:,2:end))/length(drifter_id)*100;
rej_rate_lag0 = num2str(rejection_rate(1));
disp(['Rejection rate using selected subsampling lag: ' rej_rate_lag0]);
figure('name', 'Percentage of rejected tests');
plot(distances, sum(results(:,3:end))/length(drifter_id)*100);
title('Percentage of rejected tests depending on subsampling lag');
xlabel('lag');
ylabel('Rejection percentage');