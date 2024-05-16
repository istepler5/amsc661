% AMSC 661 Final Exam Problem 1
% Consider the van der Pol oscillator x'' - mu(1-x^2)x' + x = 0. Rewrite it
% as a system of two ODEs for x and y = dx/dt. As the parameter mu
% increases, the ODE system becomes more and more stiff. Use the DIRK2
% method for integrating the system. In all cases, use (x=2,y=0) as the
% initial condition. 
%% part (a)
% Set mu = 10^2, T_max = 200, and choose time step h = 10^-3. Implement the
% DIRK2 method and integrate the ODE. Measure the CPU time. Plot x and y
% versus t in Figure 1 and y versus x in Figure 2. 
% Try to increase mu to 10^3 while leaving all other settings unchanged and
% plot the same graphs. It should be clear that the numerical solution with
% a constant time step that is small enough so that the figure in the
% xy-plane looks nice and for a time long enough so that the whole period
% is displayed will take too long.

% Note, the RHS of the ODE system is autonomous (does not depend on time).

h = 1e-3;
Tmax = 200;
t = 0:h:Tmax;
mu = 100;
g = 1 - (1/sqrt(2));    % value of gamma for DIRK2 method
init_cond = [2; 0];
tol = 1e-8;

u = zeros(2,length(t));      % keeps track of the solution over time
u(:,1) = init_cond;     % setting the initial condition

% start the time measurement of the computation
tic
% for loop to calculate u_n from u_n-1
for i = 2:length(t)
    [~,u(:,i),~] = DIRK2step(u(:,i-1),h,tol,20,mu);
end

% end the time measurement
part_a_time = toc

% plots
figure(1);

hold on;
plot(t,u(1,:));
plot(t,u(2,:));
legend('x','y');
xlabel('time'); ylabel('numerical solution');

figure(2);
plot(u(1,:),u(2,:));
xlabel('x'); ylabel("y=x'");

%% part (b)
% Make the time step adaptive using the approach illustrated on page 22 in
% ODEsolvers.pdf. Set atol = rtol = 1e-5, propose one more vector b in the
% Butcher array that makes the method 1st order, introduce the automatic 
% error estimate e and the step acceptance criterion
% ||e|| < atol + rtol · ||[x, y]||. Propose an algorithm that increases, 
% leaves the same, or decreases the time step depending on the relationship
% between ||e|| and atol + rtol · ||[x, y]||. Use the resulting solver to 
% integrate the ODE with µ = 106 on the time interval [0, 2 · 10^6]. 
% Measure the CPU time. Plot x and y versus T and y versus x.

mu = 1e6;
t0 = 0;
tmax = 2*mu;
h = 1;       % timestep to start, then will be adaptive from there
g = 1 - (1/sqrt(2));
bhat = [1-g g];
b = [0.5 0.5];
e = b - bhat;
tol = 1e-8;
atol = 1e-5;
rtol = atol;

% don't know the timesteps so don't know the sizes of the t and u vectors -
% initialize as single values and add to them with each loop
t_vec = t0;
u = [2; 0];         % initialize u with the initial condition
u_arr = u;

tic
while t_vec(end) < tmax
    % [~,u(:,end+1),k] = DIRK2step(u(:,end),h,tol,20,mu);
    [~,unew,k] = DIRK2step(u,h,tol,20,mu);
    % e_n will do e1 * k1 (column vec) and e2 * k2 (column vec) and sum
    e_n = h * norm(e(1)*k(:,1) + e(2)*k(:,2));
    eps = atol + rtol * norm(unew);

    while e_n > eps
        h = h * 0.9 * sqrt(eps/e_n);
        % redo the u from the same timestep with new h
        [~, unew, k] = DIRK2step(u,h,tol,20,mu);
        e_n = h * norm(e(1)*k(:,1) + e(2)*k(:,2));
        eps = atol + rtol * norm(unew);
        if e_n <= eps
            break
        end
    end
    t_vec(end+1) = t_vec(end) + h;
    u = unew;
    u_arr(:,end+1) = unew;
    
    if e_n < eps
        % set h for the next timestep
        h = h * 0.9 * sqrt(eps/e_n);
    end
     
end

part_b_time = toc

% plots
figure(1);
hold on;
plot(t_vec,u_arr(1,:));
plot(t_vec,u_arr(2,:));
hold off;
legend('x','y');
xlabel('time'); ylabel('numerical solution');

figure(2);
plot(u_arr(1,:),u_arr(2,:));
xlabel('x'); ylabel("y=x'");

%% functions needed

function output = f(u,mu)
% gives the RHS of the ODE system
    output = zeros(2,1);
    x = u(1);
    y = u(2);
    output(1) = y;
    output(2) = mu*(1-x^2)*y - x;
end

function output = fjac(u,mu)
% gives the Jacobian of the RHS of the ODE system
    output = zeros(2);
    x = u(1);
    y = u(2);
    output(1,1) = 0;
    output(1,2) = 1;
    output(2,1) = -2*mu*x*y - 1;
    output(2,2) = mu*(1-x^2);
end

% DIRK2 (from Homework 3)
function knew = NewtonIterDIRK2(y,h,k,gamma,mu)
    aux = y + h*gamma*k;
    F = k - f(aux,mu);
    DF = eye(2) - h*gamma*fjac(aux,mu);
    knew = k - DF\F;
end

function [ynew,y2new,k] = DIRK2step(y,h,tol,itermax,mu)
    gamma = 1.0 - (1.0/sqrt(2));
    k1 = f(y,mu);
    for j = 1 : itermax
        k1 = NewtonIterDIRK2(y,h,k1,gamma,mu);
        if norm(k1 - f(y + h*gamma*k1, mu)) < tol
            break
        end
    end
    
    % calculate first part of y2 (second-order solution) 
    % and y (first-order solution)
    y2 = y + h*(1-gamma)*k1;
    y = y + h*0.5*k1;      
    
    k2 = k1;
    for j = 1 : itermax
        k2 = NewtonIterDIRK2(y2,h,k2,gamma,mu);
        aux = y2 + h*gamma*k2;
        if norm(k2 - f(aux,mu)) < tol
            break
        end
    end
    y2new = aux;
    ynew = y + h*0.5*k2;
    k = [k1 k2];
end
