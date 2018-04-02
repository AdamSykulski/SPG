classdef AR1C
    %Class that defines a complex-valued first-order model, with fixed
    %damping and variance parameters a and sigma^2. The noise is
    %complex-valued Gaussian proper.
    
    properties
        a = 0.9;        %Damping parameter
        sigma = 1;      %Amplitude of the noise
        theta = 0.02;   %Frequency parameter when the latter is constant
        theta_t;        %Changing frequency parameter sequence
    end
    
    methods(Static)
        function c = autocov_seq(a, sigma, maxLag, includeNegative)
            if includeNegative == false
                lags = (0:maxLag);
                c = sigma^2*a.^(lags)/(1-a^2);
            else
                lags = (-maxLag:maxLag);
                c = sigma^2*a.^(abs(lags))/(1-a^2);
            end
        end
    end
    
    methods
        %Constructor
        function obj = AR1C(a, sigma, theta)
            if nargin > 0
                obj.a = a;
                obj.sigma = sigma;
                obj.theta = theta;
            end
        end
        
        %Set methods
        function obj = set.a(obj, v)
            obj.a = v;
        end
        
        function obj = set.sigma(obj, v)
            obj.sigma = v;
        end
        
        function obj = set.theta_t(obj, array)
            %Function description:
            %Used to define the sequence of frequency parameters.
            obj.theta_t = array;
        end
        
        %Simulation method.
        %Parameters:
        %   -T: sample size
        %   -nbSamples: number of independent samples
        %Returned value:
        %   -data: [TxnbSamples] double array. The simulated data with each
        %   sample in a separate column.
        function data = simulate(obj, T, nbSamples)
            data = zeros(T, nbSamples);
            eps = obj.sigma*(randn(T, nbSamples) + 1i*randn(T, nbSamples));
            data(1, :) = 1/sqrt(2*(1-obj.a^2))*eps(1, :);
            for t=2:T
                data(t,:) = obj.a*exp(1i*obj.theta_t(t))*data(t-1, :) + 1/sqrt(2)*eps(t, :);
            end
        end
    end
    
end

