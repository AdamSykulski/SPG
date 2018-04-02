for ii = 1:100 % 100 repeats
for slab_damp = [1 2 3 4 5 6 7 8]; % range of damping parameters
for jj = 1:10; % range of meridional mean flows
    [ii slab_damp jj]
    meanV = -(jj-1)/10; % range of meridional mean flows in m/s
    %% GENERATING THE DATA
    % set the initial latitude of the simulation (the longitude doesn't matter dynamically)
    initialLatitude = 35;
    % set the zonal mean flow
    meanU = 0;
    %  Number of days to compute the simulation
    totalDays = 60;
    % Resolution of the model output, in days (2/24 = 2 hours).
    outputResolution = 2/24;
    % We interpret the wind stress as a uniform force throughout the water body, 
    % so we therefore need to specify the depth of this 'slab' of water.
    slab_depth = 50;
    % Our wind data will be driven by white noise forcing, but the
    % magnitude of the forcing is determined from real data as follows
    % First load the winds from the PAPA mooring in the Gulf of Alaska
    wind = load('wind');
    % Limit the length of the wind time series to the total number of days.
    endIndex = find( wind.t > totalDays*86400, 1, 'first');
    wind.t = wind.t(1:endIndex);
    wind.u = wind.u(1:endIndex);
    wind.v = wind.v(1:endIndex);
    rho_water = 1025; % units of kg/m^3
    [tau_t, tau_x, tau_y] = StressFromWindVector( wind.t, wind.u, wind.v );
    tau_x = std(tau_x)*randn(1,length(wind.t)); % now replace with white noise, keeping same magnitude
    tau_y = std(tau_y)*randn(1,length(wind.t)); % now replace with white noise, keeping same magnitude
    % Convert the stress to a body force
    f_t = wind.t;
    f_x = tau_x / (rho_water * slab_depth);
    f_y = tau_y / (rho_water * slab_depth);
    % Convert the slab damping parameter 
    r = 1/(slab_damp*86400);
    % Initial conditions (longitude, latitude, u, v) in degrees and meters/second
    initialConditions = [0, initialLatitude, 0.0, 0.0];
    % time points at which we want the output
    outputTimes = (0:outputResolution:totalDays)*86400;
    % Now solve the differential equation
    func = @(time, inVec) ForcedInertialOscillationsOnTheEarthFlux(time, inVec, r, f_t, f_x, f_y, 0, meanU, meanV);
    options = odeset('RelTol',1e-3,'AbsTol',1e-5);
    [T, X] = ode45(func,outputTimes,initialConditions, options);
    t = T;
    tDays = t/86400;
    longitude = X(:,1);
    latitude = X(:,2);
    u_inertial = X(:,3);
    v_inertial = X(:,4);
    u_total = u_inertial + meanU;
    v_total = v_inertial + meanV;
    cv = u_inertial + sqrt(-1)*v_inertial; % this is the time series from the numerical model
    %% Matern OU fit
    Delta = 2; % TIME BETWEEN OBSERVATIONS
    ZEROF=1; % Include zero-frequency? 1=yes, 0=no;
    options=optimset('GradObj','on','MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-10,'TolX',1e-10,'Display','on'); % Fminsearch choices
    CF = -4*7.2921*10^(-5)*3600*sin(mean(latitude)/90*(pi/2)); % average Coriolis frequency
    phi_t = -4*7.2921*10^(-5)*3600*sin(latitude/90*(pi/2)); % Coriolis frequencies.
    UPCNT=1; % Fraction of positive frequencies included in fit (from 0 to 1)
    LPCNT=1; % Fraction of negative frequencies included in fit (from 0 to 1)
    N = length(cv); omega=0:2*pi/N:2*pi*(1-1/N); omega=fftshift(omega); omega(1:floor(N/2))=omega(1:floor(N/2))-2*pi; % Fourier frequencies
    MF = floor(N/2)+1; LB = round((MF-1)*(1-LPCNT)+1); UB = round(MF+UPCNT*(N-MF)); % Frequency indices in our estimation range
    SZ=(Delta/N)*(abs(fft(cv))).^2; SZ=fftshift(SZ); % Power Spectrum
    %% INITIAL PARAMETERS (same as in Sykulski et. al (2016) JRSSc for stationary drifters)
    xb=zeros(1,3);
    IF=round(N*0.5*(1+0.5*CF/pi)); 
    if CF > 0 % find location of inertial peak
        [dum1,fmax] = max(SZ(IF:N)); fmax=fmax+IF-1;
    else
        [dum1,fmax] = max(SZ(1:IF));
    end 
    NF=floor((abs(MF-fmax))/3); % Number of frequencies we now analyze over
    valmax = SZ(fmax); xb(2)=CF; % corresponding frequency peak in radians
    xb(3)=quantile(SZ([fmax-NF:fmax-1 fmax+1:fmax+NF]).*(omega([fmax-NF:fmax-1 fmax+1:fmax+NF])-xb(2))'.^2./(valmax-SZ([fmax-NF:fmax-1 fmax+1:fmax+NF])),0.5); % solve simultaneous equations for xb(1) and xb(3), take median
    xb(1)=valmax*xb(3); xb(1) = abs(xb(1))^0.5; xb(3) = abs(xb(3))^0.5; % square root of these to arrive at starting parameters
    %% LIKELIHOOD----------------------------------------------
    %%The stationary likelihood of Sykulski et al. (2016) JRSSc
    minbound = [0  1  pi/(N*xb(3))]; % minimum parameter range
    maxbound = [inf 1 inf]; % maximum parameter range
    [x1b,fval1,exitflag1]=fminsearchbnd(@(x) OUmodel(x,xb,SZ',N,LB,UB,MF,ZEROF), ones(1,3),minbound,maxbound, options); 
    x1 = x1b .* xb;
    %%Our nonstationary likelihood
    xb2 = xb; % same starting values as for stationary method
    minbound = [0  0  pi/(N*xb2(3))]; % same minimum parameter range
    maxbound = [inf 0 inf]; % same maximum parameter range
    [x2, fval2, exitflag2, time2, ker] = fitOU( SZ', phi_t, N, LB, UB, MF, ZEROF, xb2, minbound, maxbound, options);
    %% Stored output
    STAT(slab_damp,jj,ii,:)=x1;
    NSTAT(slab_damp,jj,ii,:)=x2;
end
end
end
%%
save('STAT.mat','STAT');
save('NSTAT.mat','NSTAT');