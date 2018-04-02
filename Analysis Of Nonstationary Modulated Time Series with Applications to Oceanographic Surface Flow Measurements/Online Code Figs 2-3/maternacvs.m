function [acv] = maternacvs(Q,N,Delta)
acv = zeros(1,N); % acvs from lag 0 to N-1 (not 1 to N)
acv(1)=Q(1)^2*beta(0.5,Q(2)-0.5)/(2*pi*abs(Q(3))^(2*Q(2)-1)); % Variance
acv(2:N)=Q(1)^2*(abs(Q(3))*Delta*(1:N-1)).^(Q(2)-.5).*besselk(Q(2)-.5,abs(Q(3))*Delta*(1:N-1))./(sqrt(pi)*2^(Q(2)-.5).*(abs(Q(3)))^(2*Q(2)-1)*gamma(Q(2)));