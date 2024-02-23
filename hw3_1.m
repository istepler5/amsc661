% AMSC Homework 3 Problem 1

% Show that the DIRK method with given Butcher array and gamma = 1 -
% 1/root2  is L-stable and 2nd-order accurate. 

% First, we need to show that the method is A-stable, which we can do by
% plotting the RAS of the method. 

nx = 200;
ny = 200;
x = linspace(-5,5,nx);
y = linspace(-5,5,ny);
[xg,yg] = meshgrid(x,y);

z = xg + 1i * yg;
gamma = 1 - (1/sqrt(2));

f = (z - 2*gamma*z + 1) ./ (1-gamma*z).^2;
absf = real(f).^2 + imag(f).^2;

%%
figure(1);
contourf(xg,yg,absf,[0 1]);
colorbar;
grid on;
xlabel('Re(z)'); ylabel('Im(z)');
title('RAS of DIRK2')



