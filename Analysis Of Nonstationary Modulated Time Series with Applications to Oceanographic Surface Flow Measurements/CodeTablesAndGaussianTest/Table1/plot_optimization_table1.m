%Author: Arthur P. Guillaumin
%Date: 13/06/2016
%File description:
%This script generates independent time series according to a complex-valued AR(1) process with a time-varying
%frequency parameter. The variation of this frequency parameter is generated according to a bounded random walk 
%around some value. A certain number of independent time series are generated according to this process.
%For each of these, estimation is carried out for different subsamples sizes.
%Two estimation methods are implemented.
%The first one is the blurred Whittle likelihood, where a constant frequency parameter is assumed.
%The second one is the frequency-domain likelihood that accounts for the changing frequency parameter.
%Summaries of the estimation procedure are reported in tables automatically
%generated at the end of the script.
%
%Running time < 1h on an average computer.

clear all;                  %Delete all
rng(8167,'twister');        %Set the seed and the random generator

%Parameters for the script
sample_size = 2^12;		%Total length of each sample
sample_sizes_list =  2.^(12);     %Subsample sizes used for estimation
nb_sample_sizes = length(sample_sizes_list);            %Number of sub-sample sizes
nb_independent_samples = 20;                          %Number of independent samples

%Parameters of the latent AR(1) process, to be estimated later.
r = 0.8;		%regression parameter
sigma = 1;		%Innovation amplitude


%The sequence of frequency parameters is generated according to a bounded
%random walk (equations (37)-(39))
betas = zeros(sample_size, nb_independent_samples);		%Sequences of frequency parameter
Delta = 1;
gamma = pi/2;
A = 1/20;
betas(1,:) = max(min(gamma + A*randn(1, nb_independent_samples), gamma + Delta), gamma - Delta);
for t=2:sample_size
    betas(t,:) = max(min(betas(t-1,:) + A*randn(1, nb_independent_samples), gamma + Delta), gamma - Delta);
end

%Generated data will be stored in Z
Z = zeros(sample_size, nb_independent_samples);
%Estimates stored in the following arrays
%The following are used for the stationary estimation method
r_estimates = zeros(nb_sample_sizes, nb_independent_samples);
sigma_estimates = zeros(nb_sample_sizes, nb_independent_samples);
estimate_times = zeros(nb_sample_sizes, nb_independent_samples);
exit_flags =  zeros(nb_sample_sizes, nb_independent_samples);
%The following are used for the nonstationary estimation method
r_estimates_2 = zeros(nb_sample_sizes, nb_independent_samples);
sigma_estimates_2 = zeros(nb_sample_sizes, nb_independent_samples);
estimate_times_2 = zeros(nb_sample_sizes, nb_independent_samples);
exit_flags_2 =  zeros(nb_sample_sizes, nb_independent_samples);

%Generation of the modulated complex-valued AR(1) process according to
%equation (22).
innovationsR = randn(sample_size, nb_independent_samples);	%Real part of innovations
innovationsI = randn(sample_size, nb_independent_samples);	%Imaginary part of innovations
innovations = innovationsR + 1i * innovationsI;				%Complex-valued innovations
Z(1, :) = sigma/sqrt(2*(1-r^2)) .* innovations(1,:);        
for t = 1:sample_size-1
    Z(t+1, :) = r*exp(1i*betas(t,:)).*Z(t, :) + sigma/sqrt(2)*innovations(t+1,:);
end

%Options for the optimization procedure
options=optimset('GradObj','on','MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-3,'TolX',1e-7,'Display','on');
%Starting parameters for the optimization procedure
rStart = 0.1;	
sigmaStart = 0.1;


%Estimations
%Loop over different sample sizes
for i_sample_size = 1:nb_sample_sizes
    size = sample_sizes_list(i_sample_size);            %Current sample size
    
    %Computation of periodograms. The fft function applied to a matrix
    %returns the ffts of each column of that matrix. Note that with
    %subsample according to the current sampe size.
    P = 1/size*abs(fft(Z(1:size,:))).^2;       
    
    %Loop over different samples
    for i_sample = 1:nb_independent_samples
		disp(['Sample size: ' num2str(size)])
        disp(['Sample nb ' num2str(i_sample)]);
		
        %Mean frequency over the sample. This is used to fit the stationary
        %model where we make the assumption of a constant frequency
        %parameter.
        phi_0 = mean(betas(1:size-1, i_sample));
        
        %True modulating sequence. See equation (23).
        g = exp(1i*cumsum([0 betas(1:size-1,i_sample)']));
        
		%Computing the kernels from the modulating sequences. We save the
		%computation times for both kernels.
        tic; ker0 = (1-(0:size-1)/size) .* exp(1i*phi_0*(0:size-1)); time_ker0 = toc;   %Stationary model
		tic; ker = kernel(g); time_ker = toc;                                           %Non-stationary model
		
        %Estimation via method not accounting for changing frequency
        tic
        [est,fval,exitflag] = fminsearch(@(x) lkh_(P(:, i_sample)', x(1), x(2), size, ker0), [rStart, sigmaStart], options);
        estimate_times(i_sample_size, i_sample) = toc + time_ker0;
        r_estimates(i_sample_size, i_sample) = est(1);
        sigma_estimates(i_sample_size, i_sample) = est(2);
        exit_flags(i_sample_size, i_sample) = exitflag;
        disp(est);
        
        %Estimation via method accounting for changing frequency
        %For this method we use the function fminseachbnd. This function
        %allows for constrained optimization. We use it as for sample size
        %128, 4 samples (out of 2000) led to a non-converging estimation.
        %Define the output function
        frames(500) = struct('cdata',[],'colormap',[]);
        output_func = @(x, optim_values, state)add_frame( P(:,i_sample), x, ker, frames, optim_values);
        options=optimset('GradObj','on','MaxFunEvals',100000,'MaxIter',...
            10000,'TolFun',1e-3,'TolX',1e-7,'Display','on',...
            'OutputFcn', output_func);
        tic
        [est2,fval2,exitflag2] = fminsearchbnd(@(x) lkh_(P(:, i_sample)', x(1), x(2), size, ker), [rStart, sigmaStart], [0 0], [1, inf], options);
        %[est2,fval2,exitflag2] = fminsearch(@(x) lkh_(P(:, i_sample)', x(1), x(2), size, ker), [rStart, sigmaStart], options);
        estimate_times_2(i_sample_size, i_sample) = toc + time_ker;
        r_estimates_2(i_sample_size, i_sample) = est2(1);
        sigma_estimates_2(i_sample_size, i_sample) = est2(2);
        exit_flags_2(i_sample_size, i_sample) = exitflag2;
        disp(est2);
    end
end


%We compute summaries for both estimation methods, to be reported in
%tables.

%Method 1.
bias_r = mean(r_estimates-r, 2);                        %Sample bias
var_r = var(r_estimates, 0, 2);                         %Sample variance
MSE_r = bias_r.^2 + var_r;                              %Sample MSE

bias_s = mean(sigma_estimates-sigma, 2);
var_s = var(sigma_estimates, 0, 2);
MSE_s = bias_s.^2 + var_s;

CPU_time = mean(estimate_times, 2);                     %Meam computation time

%Method 2.
bias_r2 = mean(r_estimates_2-r, 2);
var_r2 = var(r_estimates_2, 0, 2);
MSE_r2 = bias_r2.^2 + var_r2;

bias_s2 = mean(sigma_estimates_2-sigma, 2);
var_s2 = var(sigma_estimates_2, 0, 2);
MSE_s2 = bias_s2.^2 + var_s2;

CPU_time2 = mean(estimate_times_2, 2);


%%The tables are generated automatically using uitable.
summary = [bias_r';var_r';MSE_r'; bias_s';var_s';MSE_s';CPU_time'];
cnames = {'128', '256', '512', '1024', '2048', '4096'};   %Column names
cFormats = {'short E', 'short E', 'short E', 'short E', 'short E', 'short E'};
rnames={'Bias(a)', 'Variance (a)', 'MSE (a)', 'Bias(\sigma)', 'Variance (\sigma)', 'MSE (\sigma)', 'CPU time (s)'}; %Row names
f = figure('Name', 'Stationary method');
handle=uitable(f, 'Data', summary, 'ColumnName',cnames,'RowName',rnames, 'ColumnFormat', cFormats);   %User Interface table
set(handle, 'Position', get(handle, 'Extent'));



summary2 = [bias_r2';var_r2';MSE_r2'; bias_s2';var_s2';MSE_s2';CPU_time2'];
cnames = {'128', '256', '512', '1024', '2048', '4096'};   %Column names
rnames={'Bias(a)', 'Variance (a)', 'MSE (a)', 'Bias(\sigma)', 'Variance (\sigma)', 'MSE (\sigma)', 'CPU time (s)'}; %Row names
f = figure('Name', 'Nonstationary method');
handle=uitable(f, 'Data', summary2, 'ColumnName',cnames,'RowName',rnames, 'ColumnFormat', cFormats);   %User Interface table
set(handle, 'Position', get(handle, 'Extent'));