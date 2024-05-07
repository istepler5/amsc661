% AMSC 661 Homework 12 Problem 4
% Consider the Burgers Equation u_t + [0.5 u^2]_x = 0 with the initial
% condition u_0(x) = 1 on [0,1] and u_0(x) = 0 otherwise. Implement
% Lax-Friedrichs, Richtmeyer, MacCormack, and Godunov and apply them to the
% problem above. Compute the numerical solution with the same time step and
% plot it at times 0, 1, 2, 3, 4, 5, and 6. Plot the exact solution as
% well. 

% initializing stuff
k = .01;        % time-step
h = 0.01;      % space-step
x = -2:h:4-h;
t = 0:k:6;
u = zeros(length(x),length(t));     % x is contained in the rows, t in cols
u_temp = zeros(length(x),1);        % temporary column vector for later

% initial condition
ind = find(x >= 0 & x <= 1);
u(ind,1) = 1;

% change number to run different methods
method = 4;

if method == 1   % Lax-Friedrichs
    for i = 2:length(t)             % time for loop
        u(1,i) = u(1,i-1) - (k/h) * ( LF_flux(u(2,i-1),u(1,i-1),h,k) ...
                 - LF_flux(u(1,i-1),u(end,i-1),h,k) );
        
        for j = 2:length(x)-1       % space for loop
            u(j,i) = u(j,i-1) - (k/h) * ( LF_flux(u(j+1,i-1),u(j,i-1),h,k) ...
                     - LF_flux(u(j,i-1),u(j-1,i-1),h,k) );
        end
        u(end,i) = u(end,i-1) - (k/h) * ( LF_flux(u(1,i-1),u(end,i-1),h,k) ...
                 - LF_flux(u(end,i-1),u(end-1,i-1),h,k) );
    end

    % implement periodic boundary conditions
    u(end+1,:) = u(1,:);
    % for plotting later
    title_str = "Numerical and exact solution for Lax-Friedrichs";
end

if method == 2  % Richtmeyer
    for i = 2:length(t)
        u_temp = zeros(length(x),1);
        u_temp(1) = 0.5 * (u(1,i-1) + u(2,i-1)) - (0.5*k/h)...
                    * ( fun(u(2,i-1)) - fun(u(1,i-1)) );
        u_temp(end) = 0.5 * (u(end,i-1) + u(1,i-1)) - (0.5*k/h)...
                      * ( fun(u(1,i-1)) - fun(u(end,i-1)) );
        u(1,i) = u(1,i-1) - (k/h) * ( fun(u_temp(1)) - fun(u_temp(end)) );
        
        for j = 2:length(x)-1
            u_temp(j) = 0.5 * (u(j,i-1) + u(j+1,i-1)) - (0.5*k/h) ...
                        * ( fun(u(j+1,i-1)) - fun(u(j,i-1)) );
            u(j,i) = u(j,i-1) - (k/h)*(fun(u_temp(j)) - fun(u_temp(j-1)));
        end
        u(end,i) = u(end,i-1) - (k/h)*(fun(u_temp(end)) - fun(u_temp(end-1)));
    end
    
    % implement boundary conditions
    u(end+1,:) = u(1,:);
    % for plotting later
    title_str = "Numerical and exact solution for Richtmeyer";
end

if method == 3  % MacCormack
    for i = 2:length(t) % time for loop
        u_temp(1) = u(1,i-1) - (k/h) * ( fun(u(2,i-1)) - fun(u(1,i-1)) );
        u_temp(end) = u(end,i-1) - (k/h)*(fun(u(1,i-1)) - fun(u(end,i-1)));
        u(1,i) = 0.5 * (u(1,i-1) + u_temp(1)) - (0.5*k/h)...
                 * ( fun(u_temp(1)) - fun(u_temp(end)) );

        for j = 2:length(x)-1
            u_temp(j) = u(j,i-1) - (k/h)*(fun(u(j+1,i-1)) - fun(u(j,i-1)));
            u(j,i) = 0.5 * (u(j,i-1) + u_temp(j)) - (0.5*k/h)...
                     * ( fun(u_temp(j)) - fun(u_temp(j-1)) );
        end
        u(end,i) = 0.5 * (u(end,i) + u_temp(end)) - (0.5*k/h)...
                   * ( fun(u_temp(end)) - fun(u_temp(end-1)) );
    end

    % implement boundary conditions
    u(end+1,:) = u(1,:);
    % for plotting later
    title_str = "Numerical and exact solution for MacCormack";
end

if method == 4 % Godunov
   for i = 2:length(t)
       u(1,i) = u(1,i-1) - (k/h) * (Godunov_flux(u(1,i-1),u(2,i-1)) ...
                - Godunov_flux(u(end,i-1),u(1,i-1)));
       for j = 2:length(x)-1
           u(j,i) = u(j,i-1) - (k/h) * (Godunov_flux(u(j,i-1),u(j+1,i-1))...
                    - Godunov_flux(u(j-1,i-1),u(j,i-1)));
       end
       u(end,i) = u(end,i-1) - (k/h)*(Godunov_flux(u(end,i-1),u(1,i-1)) ...
                  - Godunov_flux(u(end-1,i-1),u(end,i-1)));
   end

   % implement boundary conditions
   u(end+1,:) = u(1,:);
   % for plotting later
    title_str = "Numerical and exact solution for Godunov";
end

% plotting for any of the methods
% make exact solution at t=0,1,2,3,4,5,6
exact_sol = zeros(length(x),7);
tt = 0:6;
for i = 1:3
    for j = 1:length(x)
        if x(j) >= 0 && x(j) <= tt(i)
            exact_sol(j,i) = x(j)/tt(i);
        elseif x(j) > tt(i) && x(j) < 1 + 0.5*tt(i)
            exact_sol(j,i) = 1;
        end
    end
end
for i = 4:length(tt)
    for j = 1:length(x)
        if x(j) >= 0 && x(j) <= sqrt(2*tt(i))
            exact_sol(j,i) = x(j)/tt(i);
        end
    end
end
exact_sol(end+1,:) = exact_sol(1,:);


t_ind = find(t==0 | t==1 | t==2 | t==3 | t==4 | t==5 | t==6);
figure(1);
for i = 1:length(tt)
    hold on;
    plot([x,x(end)+h], u(:,t_ind(i)), 'b-');
    plot([x,x(end)+h], exact_sol(:,i), 'r');
end
title(title_str);
xlabel('x'); ylabel('u');
qw{1} = plot(nan,'b');
qw{2} = plot(nan,'r');
legend([qw{:}],{'Numerical','Exact'});


                


%% useful functions

function output = fun(u)
    output = 0.5 * u^2;
end

function output = LF_flux(u_p,u,h,k)
    output = -(0.5*h/k) * (u_p - u) + 0.5 * (fun(u_p) + fun(u));
end

function output = Godunov_flux(uL,uR)
    if uL <= uR
        x = fminbnd(@(u) fun(u), uL, uR);
        output = fun(x);
    else
        x = fminbnd(@(u) -fun(u), uR, uL);
        output = fun(x);
    end
end
