% AMSC 661 Homework 1, Question 5(c) 

% Apply the linear two-step explicit method to the 2D gravity problem with
% a unit-circle solution. Integrate the numerical solution for two periods
% (t=0 to t=4pi). Show that the solution blows up by reporting the norm of
% the solution at time 4pi for h = 2pi/N with N = 20, 40, 80. Also, plot
% the x and y components of the solution in the xy-plane. 

% Note: since x, y, u, and v are used in the definition of the problem, I
% will use z to denote the numerical solution to prevent confusion. z_next
% will be the solution at the n+1-st timestep, z_curr is the solution at
% the current n-th timestep, and z_prev is the solution at the n-1-st
% timestep. 

N_vec = [20,40,80];
z0 = [1;0;0;1];

[norm20,x20,y20] = lin_2step_expl(20,z0);
[norm40,x40,y40] = lin_2step_expl(40,z0);
[norm80,x80,y80] = lin_2step_expl(80,z0);

figure(1);
plot(x20,y20,'r',x40,y40,'b',x80,y80,'g');
xlabel('x axis'); ylabel('y axis');
title('Solutions for linear 2-step explicit method');
legend('N=20','N=40','N=80');


function [norm_sol,x,y] = lin_2step_expl(N, z0)
% the linear two-step explicit method for the 2-D Gravity problem 
% calculated for two periods (time = 0 to 4pi)
    h = 2*pi / N;       % size of timestep
    Nt = 2*N;           % number of timesteps
    x = zeros(Nt,1);
    y = zeros(Nt,1);

    z_prev = z0;
    x(1) = z_prev(1);
    y(1) = z_prev(2);
    
    % Calculate the solution at the first timestep before entering the loop
    % done by setting z_n-1 = 0 (not sure this is the right way to do this)
    z_curr = -4*z0 + h*(4*gravity_RHS(z0));
    % z_curr = z0;
    x(2) = z_curr(1);
    y(2) = z_curr(2);

    for n = 3:Nt
        z_next = -4*z_curr + 5*z_prev + h*(4*gravity_RHS(z_curr) + ...
            2*gravity_RHS(z_prev));

        % track the x and y components of the solution for each iteration
        x(n) = z_next(1);
        y(n) = z_next(2);
        
        % shift the solutions for the next timestep
        z_prev = z_curr;
        z_curr = z_next;
    end
    
    norm_sol = norm(z_curr);

    fprintf('The norm at time 4pi for N = %d is %d.\n', N, norm_sol);
end

function f_z = gravity_RHS(z_input)
% Calculate the right hand side of the 2D gravity ODE, which does not
% depend on time (only the variables x,y,u,v). 
    x = z_input(1);
    y = z_input(2);
    u = z_input(3);
    v = z_input(4);
    f_z = [u; v; -x/(x^2+y^2); -y/(x^2+y^2)];
end
