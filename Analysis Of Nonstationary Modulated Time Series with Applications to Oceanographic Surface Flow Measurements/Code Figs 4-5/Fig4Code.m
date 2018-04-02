% The damping parameter expressed as a decay time in days.
slab_damp = 4;
for jj = 1:10;
    meanV = -(jj-1)/10; % range of meridional mean flows in m/s
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
    options = odeset('RelTol',1e-3,'AbsTol',[1e-5]);
    [T, X] = ode45(func,outputTimes,initialConditions, options);
    t = T;
    tDays = t/86400;
    longitude = X(:,1);
    latitude = X(:,2);
    u_inertial = X(:,3);
    v_inertial = X(:,4);
    AAA(:,jj)=longitude; % save down longitudes
    BBB(:,jj)=latitude; % save down latitudes
end
%% FIG 4 IN PAPER
CCC=AAA+cumsum(ones(721,10),2)/50-40.04;
scrsz = get(0,'ScreenSize'); % The Figure size
Arthur4=figure('Position',[1 1 scrsz(3)/1 scrsz(3)/2.6]);
plot(CCC(:,2:10)+1i.*BBB(:,2:10))
xlabel('longitude'); ylabel('latitude')
exportfig(Arthur4, 'Arthur4.eps', 'width', 16, 'color', 'cmyk','Fontmode','fixed','FontSize', 14); % For exporting fig into paper