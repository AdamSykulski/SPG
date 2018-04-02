%This script does the fit of the OU and Matern model to the GDP drifters.
Delta = 1/12; % TIME BETWEEN OBSERVATIONS IN HOURS
%Optimization parameters
options=optimset('GradObj','on','MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-10,'TolX',1e-10,'Display','on'); % Fminsearch choices
ZEROF=0; % Include zero-frequency? 1=yes, 0=no;
jj=1; FFF = ['a','b','c','d']; % Figure sub-headings
%Loading the data
load blurreddrifters.mat
load drifterulysses.mat
scrsz = get(0,'ScreenSize'); % The Figure size
Arthur2=figure('Position',[1 1 scrsz(3)/2 scrsz(3)/2.6])
for drifter_id =  [85 123 149 201] % these are the drifter IDs in the Figure example
    if drifter_id == 201
        X = drifterulysses.cv(1:852);
        CF = -4*7.2921*10^(-5)*3600*sin(mean(drifterulysses.lat(1:852))/90*(pi/2)); % average Coriolis frequency   
        phi_t = -4*7.2921*10^(-5)*3600*sin(drifterulysses.lat(1:852)/90*(pi/2)); % Coriolis frequencies.
        LPCNT=1.5/6; % Fraction of negative frequencies included in fit (from 0 to 1)
        UPCNT=1.5/6; 
    else
        X = blurreddrifters.cv{drifter_id};                 % Velocities time series
        CF = -4*7.2921*10^(-5)*3600*sin(mean(blurreddrifters.lat{drifter_id})/90*(pi/2)); % average Coriolis frequency  
        phi_t = -4*7.2921*10^(-5)*3600*sin(blurreddrifters.lat{drifter_id}/90*(pi/2)); % Coriolis frequencies.
        LPCNT=0.8/6; % Fraction of negative frequencies included in fit (from 0 to 1)
        UPCNT=0.8/6; % Fraction of positive frequencies included in fit (from 0 to 1)
    end
    if CF>0 % Only fit to one side (with the inertial oscillation)
        LPCNT=0;
    else
        UPCNT=0;
    end
    %Quick analysis of the sequence of Coriolis frequencies
    CFmin = min(phi_t);
    CFmax = max(phi_t);
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
    disp('Stationary likelihood estimation...');
    tic;
    minbound = [0  1  pi*sqrt(3)/(N*xb(3)) 0 0.5/xb(5) pi*sqrt(3)/(N*xb(6))]; % minimum parameter range
    maxbound = [inf 1 inf inf 2.5/xb(5) inf]; % maximum parameter range
    [x1b,fval1,exitflag1]=fminsearchbnd(@(x) maternOUmodel(x,xb,SZ',N,LB,UB,MF,ZEROF), ones(1,6),minbound,maxbound, options); 
    time1=toc;
    %%Our nonstationary likelihood
    disp('Nonstationary likelihood method estimation...');
    xb2 = xb; % same starting values as for stationary method
    minbound = [0  0  pi*sqrt(3)/(N*xb2(3)) 0 0.5/xb2(5) pi*sqrt(3)/(N*xb2(6))]; % same minimum parameter range
    maxbound = [inf 0 inf inf 2.5/xb(5) inf]; % same maximum parameter range
    [x2, fval2, exitflag2, time2, ker] = fitOUandMatern( SZ', phi_t, N, LB, UB, MF, ZEROF, xb2, minbound, maxbound, options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DISPLAY SOME OUTPUT
    x1 = x1b .* xb;
    disp('Likelihood 1 (stationary likelihood method):');
    disp(['Estimate: ', num2str(x1)]);
    disp(['CPU time: ', num2str(time1)]);
    disp('Likelihood 2 (nonstationary likelihood method):');
    disp(['Estimate: ', num2str(x2)]);
    disp(['CPU time: ', num2str(time2)]);
    %% FIGURE 2 IN PAPER
    %Colours
    c3 = [0 0 0];
    c2 = [0.5 0.5 0.5];
    c1 = [0.5 0.5 0.5];
%     c1 = 'b';
%     c2 = 'k';
%     c3 = 'r';
    subplot(2,2,jj);
    plot(omega*6/pi,10*(log10(SZ)),'color', c1); xlim([-pi pi]); % figure of periodogram
    acv=maternacvs(x1(4:6),N,1)+complexouacvs(x1(1:3),N,1); % predicted autocovariance sequence
    ESF2=2*fft(acv.*(1-(0:N-1)/N))-acv(1); ESF3=abs(real(fftshift(ESF2))); % fft of acvs with triangle kernel to find expected blurred spectrum
    ESF_=S_(ker, x2); %Nonstationary acv
    hold on; plot(omega*6/pi,10*log10(ESF3),'color',c2,'linestyle','--','linewidth',2); % stationary expected periodogram
    hold on; plot(omega*6/pi,10*log10(ESF_),'color',c3,'linestyle','--','linewidth',2); % nonstationary expected periodogram
    hold on; plot(omega(LB:UB)*6/pi,10*log10(ESF3(LB:UB)),'color', c2,'linewidth',2); % stationary expected periodogram (modelled range)
    hold on; plot(omega(LB:UB)*6/pi,10*log10(ESF_(LB:UB)),'color', c3,'linewidth',1.5); % nonstationary expected periodogram (modelled range)
    hold on; line([CF*6/pi CF*6/pi],[-10 60],'color','k'); % average Coriolis
    hold on; line([CFmin*6/pi CFmin*6/pi] , [-10 60],'color',c2,'linestyle','--'); % minimum Coriolis
    hold on; line([CFmax*6/pi CFmax*6/pi] , [-10 60],'color',c2,'linestyle','--'); % maximum Coriolis
    xlabel('\omega (cpd)', 'FontSize', 19); ylabel('PSD (10 log_{10} m^2s^{-2} cpd^{-1})', 'FontSize', 19);
    set(findall(gcf,'type','text'),'fontSize',19,'fontWeight','normal');
    if drifter_id == 201
        xlim([-1.5 1.5]);
    else
        xlim([-1 1]); 
    end
    ylim([-10 60]);
    title(['(' num2str(FFF(jj)) ')']);  jj=jj+1;
end
%%
exportfig(Arthur2, 'Arthur2.eps', 'width', 16, 'color', 'cmyk','Fontmode','fixed','FontSize', 19); % For exporting fig into paper