%Author: Arthur P. Guillaumin
%Date: 14/06/2016
%File description:
%Simulation of a stationary real-valued AR(1) process with given damping and noise
%amplitude parameters, subject to missing values.
%For a given time point, the probability of a missing value is p, where p 
%is generated according to a cosine function. Missing values occur independently. 
%
%Simulation and estimation are carried out for different sample sizes. The
%estimation method is described in section 5 of the Journal of Time Series
%Analysis paper with doi 10.1111/jtsa.12244.
%We report summaries of the estimations in a table generated automatically.
%
%Running time < 1hour on an average computer.
clear all;
rng(212,'twister');
nbSamples = 2000;                                       %Nb of independent samples for each sample size
sampleSizes = [128 512 1024 2048 4096 8192 16384];      %Different sample sizes 
nbSampleSizes = length(sampleSizes);

%Parameters for the latent AR(1) process
a = 0.8;                                                %damping parameter
sigma = 1;                                              %noise variance

%Options for the minimization
options=optimset('GradObj','on','MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-5,'TolX',1e-9,'Display','on');

%%Inference-----------------------------------------------------
%Table in which we will store estimates and computation times
est_ = zeros(nbSampleSizes, nbSamples, 3);

%Loop over different sample sizes
for iSampleSize=(1:nbSampleSizes)
    N = sampleSizes(iSampleSize);                       %Length of the time series

    %%Data from latent process
    X = AR1R(a, sigma, N, nbSamples);
    X = X';
    
    %Corresponding modulating sequences.
    g = PeriodicBernoulliMissing(N, nbSamples, 0.5, 0.5, 10);

    %%Modulated process
    Y = g .* X;
%     for iSample=1:nbSamples
%         Y(iSample,:) = g(:, iSample)' .* X(iSample,:);
%     end

    %Periodograms
    Q = 1/N*abs(fft(Y)).^2;
    
    %Loop over different samples
    for iSample=1:nbSamples
        %Outputs the current sample being processed
        disp(['->Sample ', num2str(iSample), ', sample Size ' num2str(sampleSizes(iSampleSize))])
        tic;    %Starts the timer
        %Computation of the kernel
        ker = kernel(g(:, iSample));
        %Estimation of the parameters of the latent process
        est_(iSampleSize, iSample, 1:2) = fminsearch(@(x) lkh_(Q(:,iSample)', x(1), x(2), N, ker), [0.7, 1.1], options);
        time = toc;     %Stops the timer
        est_(iSampleSize, iSample,  3) = time;
        disp(num2str(est_(iSampleSize, iSample, 1:3)));
    end
end


%%Summary of estimation data
%r estimate
bias_r = mean(est_(:,:,1)')-a;                          %Sample bias
var_r = var(est_(:,:,1)');                              %Sample variance
MSE_r = bias_r.^2 + var_r;                              %Sample MSE

%sigma estimate
bias_s = mean(est_(:,:,2)')-sigma;
var_s = var(est_(:,:,2)');
MSE_s = bias_s.^2 + var_s;

%CPU time
CPU_time = mean(est_(:,:,3)');

%%The table is generated automatically using uitable.
summary = [bias_r;var_r;MSE_r; bias_s;var_s;MSE_s;CPU_time];
cnames = {'128', '512', '1024', '2048', '4096', '8192', '16384'};   %Column names
rnames={'Bias(a)', 'Variance (a)', 'MSE (a)', 'Bias(\sigma)', 'Variance (\sigma)', 'MSE (\sigma)', 'CPU time (s)'}; %Row names
cFormats = {'short E', 'short E', 'short E', 'short E', 'short E', 'short E', 'short E'};
f = figure;
handle = uitable(f, 'Data', summary, 'ColumnName',cnames,'RowName',rnames, 'ColumnFormat', cFormats);   %User Interface table
set(handle, 'Position', get(handle, 'Extent'))

