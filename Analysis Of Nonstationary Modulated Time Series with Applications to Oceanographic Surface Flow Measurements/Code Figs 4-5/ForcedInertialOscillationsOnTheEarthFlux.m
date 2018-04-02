% current time, in seconds
% the inVector includes (longitude,latitude,u,v) in units of degrees and meters per second
% force_t is the time series (in seconds) corresponding to the force values force_x, force_y in units of m/s^2
% r is the damping parameter in units of inverse seconds.
function outVector = ForcedInertialOscillationsOnTheEarthFlux( time, inVector,  r, force_t, force_x, force_y, mean_t, mean_u, mean_v )
	% problem constants
	R = 6384000;
	Omega = 7.292e-5;
	
	% linearly interpolate the provided forcing to the current time.
	if (length(force_t) > 1)
		f_x = interp1( force_t, force_x, time );
		f_y = interp1( force_t, force_y, time );
	else
		f_x = force_x;
		f_y = force_y;
	end
	
	% linearly interpolate the provided background flow to the current time.
	if (length(mean_t) > 1)
		meanU = interp1( mean_t, mean_u, time );
		meanV = interp1( mean_t, mean_v, time );
	else
		meanU = mean_u;
		meanV = mean_v;
	end
	
	% write out the input vector for clarity.
	longitude = inVector(1);
	latitude = inVector(2);
	u = inVector(3);
	v = inVector(4);
	
	% precompute values that will be used more than once
	f = 2 * Omega * sind(latitude);
	gamma2 = u*tand(latitude)/R;
	rad2deg = 180/pi;
	
	% now compute the flux
	outVector = zeros(4,1);
	outVector(1) = rad2deg*(u+meanU)/(R*cosd(latitude));
	outVector(2) = rad2deg*(v+meanV)/R;
	outVector(3) = (f + gamma2)*v + f_x - r*u;
	outVector(4) = -(f + gamma2)*u + f_y - r*v;
	
	
% 	f = 2 * Omega * sind(latitude);
% 	gamma2 = (u+meanU)*tand(latitude)/R;
% 	rad2deg = 180/pi;
% 	
% 	outVector = zeros(4,1);
% 	outVector(1) = rad2deg*(u+meanU)/(R*cosd(latitude));
% 	outVector(2) = rad2deg*(v+meanV)/R;
% 	outVector(3) = (f + gamma2)*(v+meanV) + f_x - r*u;
% 	outVector(4) = -(f + gamma2)*(u+meanU) + f_y - r*v;
end