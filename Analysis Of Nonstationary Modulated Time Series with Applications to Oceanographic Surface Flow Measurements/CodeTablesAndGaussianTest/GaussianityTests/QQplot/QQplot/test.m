N = 1000;
X = randn(N,1);
P = 1/N*abs(fft(X)).^2;

S = ones(N,1);

V = P./S;
U = 1-exp(-V);
U_sorted = sort(U);
figure
plot((0:N-1)/(N-1), U_sorted);