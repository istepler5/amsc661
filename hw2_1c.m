% AMSC 661 Homework 2 Problem 1(c)

% Test appropriate ODE solvers on the Arenstorf problem. Parameters and
% initial condition values come from the note in the homework document.
% First set the t_max equal to the period, compute the periodic orbit using
% DOPRI5(4) and plot it (x vs. y). Then set T_max = 100. Set error
% tolerance atol = rtol = 10^-12. Compute the numerical solution using
% ode45, ode78, and ode89. Measure CPU times and plot the orbit of the
% satellite for each solver. Comment on the CPU times and plots that you
% observe and compare the solvers. 

y0 = [0.994; 0; 0; -2.00158510637908252240537862224];
Tmax_pd = 17.0652165601579625588917206249;

%% First, compute the periodic orbit using ode45 and plot it. 
[~,results] = ode45(@arenstorf, [0;Tmax_pd],y0);

figure(1);
plot(results(:,1),results(:,2));
xlabel('x'); ylabel('y');
title('Periodic Orbit for Arenstorf Problem Using ode45');

%% Compute orbits for T_max = 100 for ode45, ode78, ode89
Tmax = 100;
error = 1e-12;
options = odeset('AbsTol',error,'RelTol',error);

t45_start = cputime;
[~,y45] = ode45(@arenstorf, [0;Tmax],y0,options);
t45_end = cputime - t45_start;
fprintf('The CPU time for ode45 is %f.\n',t45_end);

figure(2);
plot(y45(:,1),y45(:,2));
xlabel('x'); ylabel('y');
title(['Arenstorf solution using ode45, Tmax = ',num2str(Tmax)]);

t78_start = cputime;
[~,y78] = ode78(@arenstorf, [0;Tmax],y0,options);
t78_end = cputime - t78_start;
fprintf('The CPU time for ode78 is %f.\n',t78_end);

figure(3);
plot(y78(:,1),y78(:,2));
xlabel('x'); ylabel('y');
title(['Arenstorf solution using ode78, Tmax = ',num2str(Tmax)]);

t89_start = cputime;
[~,y89] = ode89(@arenstorf, [0;Tmax],y0,options);
t89_end = cputime - t89_start;
fprintf('The CPU time for ode89 is %f.\n',t89_end);

figure(4);
plot(y89(:,1),y89(:,2));
xlabel('x'); ylabel('y');
title(['Arenstorf solution using ode89, Tmax = ',num2str(Tmax)]);

%% ode function
function dydt = arenstorf(t,y)
    dydt = zeros(4,1);
    mu = 0.012277471;       % mass of the Moon
    mu1 = 1 - mu;           % mass of the Earth

    % I'm using the system from the PC9: Hamiltonian systems and
    % symplectic integrators PDF document, except r1 and r2 to 3/2 power
    % y(1) = x, y(2) = y, y(3) = x', y(4) = y'
    r1 = ((y(1) + mu)^2 + y(2)^2) ^ 1.5;
    r2 = ((y(1) - mu1)^2 + y(2)^2) ^ 1.5;

    dydt(1) = y(3);
    dydt(2) = y(4);
    dydt(3) = y(1) + 2*y(4) - (mu1*(y(1) + mu)/r1) - (mu*(y(1) - mu1)/r2);
    dydt(4) = y(2) - 2*y(3) - (mu1 * y(2) / r1) - (mu * y(2) / r2);
end