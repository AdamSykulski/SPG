function [] = plot_empirical_cdf(Y, title_fig)
%This function plots the empirical distribution of the passed Y.
 [f,x_values] = ecdf(Y);
 F = plot(x_values,f);
 set(F,'LineWidth',2);
 hold on;
 G = plot(x_values,normcdf(x_values,0,1),'r-');
 set(G,'LineWidth',2);
 title(title_fig);
 hold on;
 kde(Y);
 legend('Empirical CDF','Standard Normal CDF', 'kde', ...
     'Location','SE');
 axis([-5 5 0 1]);
end

