%Author: Arthur P. Guillaumin
%Date: 15/06/2016
%File description: 
%We generate a complex-valued AR(1) process with a spin parameter that
%varies linearly with time. 5000 independent samples of length 512 are generated.
%The estimation concerns the parameter of the AR(1) model as well as the
%parameters of the linear variation (2 parameters: mean and spread) of the
%frequency. 
%Estimation is carried out via three methods:
%-An exact likelihood, for the sake of comparison. This is easy to implement
%here due to the Markov property of the process.
%-The de-biased whittle likelihood, which assumes a constant frequency
%parameter. This parameter is estimated. Note that the spread of the
%frequency parameter cannot be estimated via this method.
%-A frequency domain pseudo-likelihood which accounts for the
%linearly-varying frequency parameter.
%
%Running time < 1/2 hour on an average computer.
rng(212, 'twister');
clear all;                                              %Delete all
nbSamples = 5000;                                       %Number of independent samples
T = 512;                                                %Length of time series

%Options for the optimization procedure
options=optimset('GradObj','on','MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-10,'TolX',1e-10,'Display','on'); % some options for fminsearch, Gradient method used

%Frequency parameter sequence. See (31) in the paper.
gamma = pi/2;                                           %Average frequency parameter
delta = 2*pi/3;                                         %Spread of frequencies
theta_t = linspace(gamma-delta/2, gamma+delta/2, T);

%%Generation of the samples
r = 0.9;
sigma = 10;
model = AR1C(r, sigma, 0.8);
model.theta_t=theta_t;
data=model.simulate(T, nbSamples);
plot(data(:,1));

%Arrays in which we save estimates. Structure will be:
%[r, sigma, gamma, delta, 
est1 = zeros(nbSamples, 7);
est2 = zeros(nbSamples, 6);                             %We don't estimate delta for the second estimation procedure.
est3 = zeros(nbSamples, 7);
startparams_save = zeros(4,nbSamples);

minBound = [0.0 0 -pi -Inf];
maxbound = [1 Inf pi Inf];

%Computation of all periodograms. 
Ps = 1/T*abs(fft(data)).^2;

%Loop over all samples
for i=1:nbSamples
    disp('---------------------------------------------');
    disp(['n. sample: ' num2str(i)]);
    
    %Periodogram of sample i
    datai = data(:,i);
    P = Ps(:,i);
    
    %Initial guesses
    [Smax, omegaMax] = max(P);
    omegaMax = median(find(P>Smax/1.5))/T*2*pi;
    omegaHalf = find(P>Smax/2);
    omegaHalf = omegaHalf(1)*2*pi/T;
    omegaQuarter = find(P>Smax/4);
    omegaQuarter = omegaQuarter(1)*2*pi/T;
    omegaMax = gamma;
    rhoI = (2-cos(omegaHalf - omegaMax)-sqrt((2-cos(omegaHalf - omegaMax))^2-1))/2;
    rhoI2 = (4-cos(omegaQuarter - omegaMax)-sqrt((4-cos(omegaQuarter - omegaMax))^2-9))/3;
    sigmaI = sqrt(Smax)*(1-rhoI);
    startParams = [rhoI sigmaI omegaMax pi/4];
    startparams_save(:,i) = startParams;
    
    %%Estimation of parameters via exact MLE
    disp('Exact likelihood...');
    tic;
    [est1(i, 1:4), est1(i, 5), est1(i, 6)] = fminsearchbnd(@(x) exactLKH(x(1), x(2), x(3), x(4), datai), startParams, minBound, maxbound ,options);
    est1(i, 7) = toc;
    disp(num2str(est1(i,:)));

    %%Frequency domain MLE (stationary)
    disp('Stationary pseudo likelihood...');
    tic;
    [est2(i, 1:3), est2(i, 4), est2(i, 5)] = fminsearchbnd(@(x) freqLKHstationary(x(1), x(2), x(3), P), [rhoI sigmaI omegaMax], [0.0 0 -pi], [1 Inf pi], options);
    est2(i, 6) = toc;
    disp(num2str(est2(i,:)));

    %%Frequency-domain MLE
    disp('Nonstationary pseudo likelihood...');
    tic;
    [est3(i, 1:4), est3(i, 5), est3(i, 6)] = fminsearchbnd(@(x) freqLKH(x(1), x(2), x(3), x(4), P), startParams, minBound,  maxbound , options);
    est3(i, 7) = toc;
    disp(num2str(est3(i,:)));
end

%Summaries of estimations
bias_1 = mean(est1(:,1:4)) - [r sigma gamma delta];
var_1 = var(est1(:,1:4));
MSE_1 = bias_1.^2 + var_1;


bias_2 = [mean(est2(:,1:3)) - [r sigma gamma] 0];
var_2 = [var(est2(:,1:3)) 0];
MSE_2 = bias_2.^2 + var_2;

bias_3 = mean(est3(:,1:4)) - [r sigma gamma delta];
var_3 = var(est3(:,1:4));
MSE_3 = bias_3.^2 + var_3;


%Table generation
table_data = [bias_1;var_1;MSE_1;bias_2;var_2;MSE_2;bias_3;var_3;MSE_3];
cnames = {'r', '\sigma', '\gamma', '\Delta'};   %Column names
cFormats = {'short E', 'short E', 'short E', 'short E'};
rnames={'Bias1', 'Variance1', 'MSE1', 'Bias2', 'Variance2', 'MSE2', 'Bias3', 'Variance3', 'MSE3'}; %Row names
f=figure;
handle=uitable(f, 'Data', table_data, 'ColumnName',cnames,'RowName',rnames, 'ColumnFormat', cFormats);   %User Interface table
set(handle, 'Position', get(handle, 'Extent'));

