%This script does the fit of the OU and Matern model.
Delta = 2; % TIME BETWEEN OBSERVATIONS
%Optimization parameters
options=optimset('GradObj','on','MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-10,'TolX',1e-10,'Display','on'); % Fminsearch choices
ZEROF=0; % Include zero-frequency? 1=yes, 0=no;
%Loading the data
load blurreddrifters.mat
for drifter_id = 1:200
    drifter_id
    X = blurreddrifters.cv{drifter_id}; % Velocities time series
    CF = -4*7.2921*10^(-5)*3600*sin(mean(blurreddrifters.lat{drifter_id})/90*(pi/2)); % average Coriolis frequency
    phi_t = -4*7.2921*10^(-5)*3600*sin(blurreddrifters.lat{drifter_id}/90*(pi/2)); % Coriolis frequencies.
    LPCNT=0.8/6; % Fraction of negative frequencies included in fit (from 0 to 1)
    UPCNT=0.8/6; % Fraction of positive frequencies included in fit (from 0 to 1)
    if CF>0 % Only fit to one side (with the inertial oscillation)
        LPCNT=0;
    else
        UPCNT=0;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRELIMINARIES
    N = length(X); omega=0:2*pi/N:2*pi*(1-1/N); omega=fftshift(omega); omega(1:floor(N/2))=omega(1:floor(N/2))-2*pi; % Fourier frequencies
    MF = floor(N/2)+1; LB = round((MF-1)*(1-LPCNT)+1); UB = round(MF+UPCNT*(N-MF)); % Frequency indices in our estimation range
    SZ=(Delta/N)*(abs(fft(X))).^2; SZ=fftshift(SZ); % Power Spectrum
    %% INITIAL PARAMETERS (same as in Sykulski et. al (2016) JRSSc for stationary drifters)
    xb=zeros(1,6);
    valMF = SZ(MF); xb(5)=1; valMAX = max(SZ); IF=round(N*0.5*(1+0.5*CF/pi)); % set slope xb(5) to 1 
    if CF > 0 % find location of inertial peak
        [dum1,fmax] = max(SZ(IF:N)); fmax=fmax+IF-1;
    else
        [dum1,fmax] = max(SZ(1:IF));
    end 
    NF=floor((abs(MF-fmax))/3); % Number of frequencies we now analyze over
    xb(6)=quantile(SZ([MF-NF:MF-1 MF+1:MF+NF]).*(omega([MF-NF:MF-1 MF+1:MF+NF]))'.^2./(valMAX-SZ([MF-NF:MF-1 MF+1:MF+NF])),0.5); % solve simultaneous equations for xb(4) and xb(6), take median
    xb(4)=valMAX*xb(6); xb(4) = abs(xb(4))^0.5; xb(6) = abs(xb(6))^0.5; % square root of these to arrive at starting parameters
    valmax = SZ(fmax); xb(2)=CF; % corresponding frequency peak in radians
    xb(3)=quantile(SZ([fmax-NF:fmax-1 fmax+1:fmax+NF]).*(omega([fmax-NF:fmax-1 fmax+1:fmax+NF])-xb(2))'.^2./(valmax-SZ([fmax-NF:fmax-1 fmax+1:fmax+NF])),0.5); % solve simultaneous equations for xb(1) and xb(3), take median
    xb(1)=valmax*xb(3); xb(1) = abs(xb(1))^0.5; xb(3) = abs(xb(3))^0.5; % square root of these to arrive at starting parameters
    %% LIKELIHOOD1----------------------------------------------
    %%The stationary likelihood of Sykulski et al. (2016) JRSSc
    minbound = [0  1  pi*sqrt(3)/(N*xb(3)) 0 0.5/xb(5) pi*sqrt(3)/(N*xb(6))]; % minimum parameter range
    maxbound = [inf 1 inf inf 2.5/xb(5) inf]; % maximum parameter range
    [x1b,fval1,exitflag1]=fminsearchbnd(@(x) maternOUmodel(x,xb,SZ',N,LB,UB,MF,ZEROF), ones(1,6),minbound,maxbound, options); 
    x1 = x1b .* xb;
    %%Our nonstationary likelihood
    xb2 = xb; % same starting values as for stationary method
    minbound = [0  0  pi*sqrt(3)/(N*xb2(3)) 0 0.5/xb2(5) pi*sqrt(3)/(N*xb2(6))]; % same minimum parameter range
    maxbound = [inf 0 inf inf 2.5/xb(5) inf]; % same maximum parameter range
    [x2, fval2, exitflag2, time2, ker] = fitOUandMatern( SZ', phi_t, N, LB, UB, MF, ZEROF, xb2, minbound, maxbound, options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stat.ml(drifter_id)=fval1;
    nstat.ml(drifter_id)=fval2;
    stat.param(drifter_id,:)=x1;
    nstat.param(drifter_id,:)=x2;
end
%% FIGURE 3 IN PAPER
scrsz = get(0,'ScreenSize'); % The Figure size
Arthur3=figure('Position',[1 1 scrsz(3)/1 scrsz(3)/2.5]); subplot(1,2,1)
scatter(1./(stat.param(:,3)*12),1./(nstat.param(:,3)*12)) % Scatter plot of damping timescales
axis square; xlim([0 10]); ylim([0 10]); xlabel('1/\lambda (days, stationary fit, median = 1.29)'); ylabel('1/\lambda (days, nonstationary fit, median = 3.42)'); 
title('(a)'); box on
hold on; line([0 10],[0 10]);
subplot(1,2,2); histogram(stat.ml-nstat.ml,50); title('(b)'); ylabel('frequency'); xlabel('l_M(\theta_M) - l_W(\theta_W)'); % Histogram of likelihood differences
exportfig(Arthur3, 'Arthur3.eps', 'width', 16, 'color', 'cmyk','Fontmode','fixed','FontSize', 16); % For exporting fig into paper
%%
disp(['Median damping timescale (stationary model): ', num2str(median(1./(stat.param(:,3)*12)))]);
disp(['Median damping timescale (nonstationary model): ', num2str(median(1./(nstat.param(:,3)*12)))]);