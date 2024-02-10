% AMSC 661 Homework 1, Question 5(d)

% Apply the midpoint rule with the Forward Euler predictor to the 2D
% gravity problem with a unit-circle solution. Use the time interval 8pi,
% set h = 2pi/N with N=20, 40, and 80 and plot the numerical solutions.

% Note: since x, y, u, and v are used in the definition of the problem, I
% will use z to denote the numerical solution to prevent confusion. z_next
% will be the solution at the n+1-st timestep, and z_curr is the solution at
% the current n-th timestep.

N_vec = [20,40,80];
z0 = [1;0;0;1];

[norm20,x20,y20] = midpt_FE_predictor(20,z0);
[norm40,x40,y40] = midpt_FE_predictor(40,z0);
[norm80,x80,y80] = midpt_FE_predictor(80,z0);

figure(1);
plot(x20,y20,'r',x40,y40,'b',x80,y80,'g');
xlabel('x axis'); ylabel('y axis');
title('Solutions for midpoint rule with Forward Euler predictor');
legend('N=20','N=40','N=80');
axis equal;


function [norm_sol,x,y] = midpt_FE_predictor(N, z0)
% the midpoint rule with Forward Euler predictor for the 2-D Gravity 
% problem calculated for 4 periods (time = 0 to 8pi)
    h = 2*pi / N;       % size of timestep
    Nt = 4*N;           % number of timesteps
    x = zeros(Nt,1);
    y = zeros(Nt,1);

    z_curr = z0;
    x(1) = z_curr(1);
    y(1) = z_curr(2);

    for n = 2:Nt
        z_next = z_curr + h * gravity_RHS(z_curr + ...
            (h/2)*gravity_RHS(z_curr));

        % track the x and y components of the solution for each iteration
        x(n) = z_next(1);
        y(n) = z_next(2);
        
        % shift the solutions for the next timestep
        z_curr = z_next;
    end
    
    norm_sol = norm(z_curr);

    fprintf('The norm at time 8pi for N = %d is %d.\n', N, norm_sol);
end

function f_z = gravity_RHS(z_input)
% Calculate the right hand side of the 2D gravity ODE, which does not
% depend on time (only on the variables x,y,u,v).
    x = z_input(1);
    y = z_input(2);
    u = z_input(3);
    v = z_input(4);
    f_z = [u; v; -x/(x^2+y^2); -y/(x^2+y^2)];
end
