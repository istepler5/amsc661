% AMSC 661 Homework 2 Problem 2

% On the same coordinate plane, plot the regions of absolute stability for
% the following ERK methods: Forward Euler, Midpoint Rule with Euler
% predictor, Kutta's method, the standard 4-stage, 4th order RK method, and
% DOPRI5(4) (use y-hat for the method of order 5).

nx = 200;
ny = 200;
x = linspace(-5,5,nx);
y = linspace(-5,5,ny);
[xg,yg] = meshgrid(x,y);

z = xg + 1i * yg;

% For each method, R(z) will be labeled f1, f2, and so on.
% We then take |R(z)| (e.g. absf1) because we will plot
% the level set |R(z)|<1 to show the RAS.

% Forward Euler
f1 = 1 + z;
absf1 = real(f1).^2 + imag(f1).^2;

% Midpoint Rule with Forward Euler predictor
f2 = 1 + z + 0.5*z.^2;
absf2 = real(f2).^2 + imag(f2).^2;

% Kutta's method
f3 = 1 + z + 0.5*z.^2 + (1/6)*z.^3;
absf3 = real(f3).^2 + imag(f3).^2;

% 4-stage, fourth order Runge-Kutta method
f4 = 1 + z + 0.5*z.^2 + (1/6)*z.^3 + (1/24)*z.^4;
absf4 = real(f4).^2 + imag(f4).^2;

% DOPRI5(4)
b = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
A = [0, 0, 0, 0, 0, 0, 0;
     0.2, 0, 0, 0, 0, 0, 0;
     3/40, 9/40, 0, 0, 0, 0, 0;
     44/45, -56/15, 32/9, 0, 0, 0, 0;
     19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0, 0;
     9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0, 0;
     35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0];
c = ones(7,1);

syms a
expr = 1 + a * b * (eye(7) + a * A + a.^2 * A^2 + a.^3 * A^3 + a.^4 * A^4 ...
       + a.^5 * A^5 + a.^6 * A^6 + a.^7 * A^7) * c;
f5_expr = simplify(expr);
f5 = double(subs(f5_expr,a,z));
absf5 = real(f5).^2 + imag(f5).^2;

%% Plotting
figure(1);

subplot(2,3,1)
contourf(xg,yg,absf1,[0 1]);
grid on;
xlabel('Re(z)'); ylabel('Im(z)');
title('Forward Euler');

subplot(2,3,2);
contourf(xg,yg,absf2,[0,1]);
grid on;
xlabel('Re(z)'); ylabel('Im(z)');
title('Midpoint Rule with FE Predictor')

subplot(2,3,3);
contourf(xg,yg,absf3,[0 1]);
grid on;
xlabel('Re(z)'); ylabel('Im(z)');
title('Kutta Method');

subplot(2,3,4);
contourf(xg,yg,absf4,[0,1]);
grid on;
xlabel('Re(z)'); ylabel('Im(z)');
title('4-stage 4th-order Runge-Kutta');

subplot(2,3,5);
contourf(xg,yg,absf5,[0,1]);
grid on;
xlabel('Re(z)'); ylabel('Im(z)');
title('DOPRI5(4)');



