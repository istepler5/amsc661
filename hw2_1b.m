% AMSC 661 Homework 2 Problem 1(b)

% Test the ODE solvers ode45 and ode15s on the Van Der Pol oscillator:
% y1' = y2, y2' = mu((1-y1^2)y2) - y1, with mu = 10, 10^2, 10^3. The
% greater mu is, the stiffer the problem. Set tmax = 1000.0. Set error
% tolerances epsilon = 10^-6, 10^-9, and 10^-12. Pick the initial condition
% y1(0) = 2, y2(0) = 0. Plot the solution on the phase plane (y1, y2). FOr
% each value of mu, make a single plot of the solution. Observe how the
% limit cycle changes as mu increases. Measure the CPU time T_CPU required
% to compute one cycle for each of the solvers for each of the values mu
% and plot log epsilon versus log CPUtime. Comment on how the CPU time
% depends on the error tolerance for each value of mu for each method. Try
% to explain your observations. 

mu = [10, 100, 1000];
index = [1,3,5];

for j = 1:length(mu)
    ode_comp(mu(j),index(j));
end

function [] = ode_comp(mu,i)
% this function will compute 3 solutions for ode45 with the three error
% values and 3 solutions for ode15s with the three error values, then
% create plots to compare the phase planes of the solutions and the errors
% versus the CPU times. 
    error = [1e-6, 1e-9, 1e-12];
    y0 = [2;0];
    tmax = 1000.0;

    for k = 1:length(error)
        options = odeset('AbsTol',error(k),'RelTol',error(k));
        t45_start = cputime;
        [~,y45] = ode45(@(t,y) vanderpol(t,y,mu),[0;tmax],y0,options);
        t45_end = cputime - t45_start;

        t15s_start = cputime;
        [~,y15s] = ode15s(@(t,y) vanderpol(t,y,mu),[0;tmax],y0,options);
        t15s_end = cputime - t15s_start;

        fprintf(['The CPU time for ode45 with mu = %d and error = %d is ' ...
            '%f.\n'],mu,error(k),t45_end);
        fprintf(['The CPU time for ode15s with mu = %d and error = %d is ' ...
            '%f.\n'],mu,error(k),t15s_end);

        figure(i);
        plot(y45(:,1),y45(:,2));
        hold on;
        plot(y15s(:,1),y15s(:,2));
        title(['Comparison of ode45 to ode15s for mu = ',num2str(mu)]);
        xlabel('y1'); ylabel('y2');
        legend('y45 error 1e-6','y15s error 1e-6','y45 error 1e-9',...
            'y15s error 1e-9','y45 error 1e-12','y15s error 1e-12',...
            Location='southeast');

        figure(i+1);
        loglog(error(k),t45_end,'o',error(k),t15s_end,'o');
        hold on;
        axis([0 1e-5 0 3]);
        title(['log(error) vs. log(CPU Time) for mu = ',num2str(mu)]);
        xlabel('log(error)'); ylabel('log(CPUtime)');
        grid on;
        legend('y45 error 1e-6','y15s error 1e-6','y45 error 1e-9',...
            'y15s error 1e-9','y45 error 1e-12','y15s error 1e-12');
    end
    hold off; hold off;
end

function dydt = vanderpol(t,y,mu)
% calculates the right hand side of the Van Der Pol oscillator system    
    dydt = zeros(2,1);
    dydt(1) = y(2);
    dydt(2) = mu * ((1-y(1)^2)*y(2)) - y(1);
end

