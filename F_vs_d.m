% displacement vs force values, copied in bc im a gremlin
d = [0, 0.006692379, 0.013384758, 0.020077137, ...
     0.026769515, 0.033461894, 0.040154273, ...
     0.046846652, 0.053539031, 0.06023141,0.066923789];

F = [0, 0.005046889, 0.016769448, ...
      0.036419621, 0.063730441, 0.099828104, ...
      0.144767553, 0.199078755, 0.263437554, ...
      0.338697924, 0.424681993];

d_fake = linspace(0,0.066923789,50);
% F=ax^2+bx least-squares fit params, according to scipy optimize
%a = 98.52325362043621;
%b = -0.3033535333676083;
%F_hat = a*d_fake.*d_fake + b*d_fake;

% F=Ax^n fit params, according to log-log lin reg in excel
n = 1.9397;
A = exp(4.3254);
F_hat = A*d_fake.^n;


plot(d_fake,F_hat);
hold on
scatter(d,F)
hold off
prettyplot('','d/a','free energy (arb.)')
