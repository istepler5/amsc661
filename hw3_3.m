% AMSC 661 Homework 3 Problem 3

% Examine the performance of DIRK methods of orders 2 and 3 on the
% Prothero-Robinson problem, plot the maximum absolute error as a function
% of the time step, and observe the two orders of convergence for each
% method, one for large hL and one for small hL. 
% Prothero-Robinson problem: y' = -L(y-phi(t)) + phi'(t), y(0) = y_0 with L
% = 10^4 and phi(t) = sin(t + pi/4). Time interval 0 to Tmax = 10. 

%% part (a)
% Pick initial condition y(0) = sin(pi/4). Compute the numerical solution
% using DIRK2 on the interval [0, Tmax] with time step h for each h from
% the following set: h = 10^-p, p in {1, 1+d, 1+2d,...6}, d = 5/24. Plot
% numerical error e(h) = max_0 to T_max {|u_n - y(t_n)|} vs. h. Use log-log
% scale. Observe error decay e = C1h for large values of h and e = C2h^2 
% for small values of h. Do the same for DIRK of order 3. What orders of 
% error decay do you observe? Also, plot reference lines. 
Tmax = 10;
L = 10^4;
p = 1:(5/24):6;
h = 10 .^ (-p);
y0_a = sin(pi/4);
gamma1 = 1 - (1/sqrt(2));
gamma2 = 0.5 + sqrt(3)/6;

% Calculate the e(h) for DIRK2 and DIRK3 for each h in the h vector
error_dirk2a = zeros(1,length(h));
error_dirk3a = zeros(1,length(h));
for j = 1:length(h)
    [error_dirk2a(j),~] = DIRK2(gamma1,Tmax,L,h(j),y0_a);
    [error_dirk3a(j),~] = DIRK3(gamma2,Tmax,L,h(j),y0_a);
end

figure(1);
c1 = 0.00004;
c2 = 0.1;
c3 = 1000;
e1 = c1*h;
e2 = c2*h.^2;
e3 = c3*h.^3;
loglog(h,e1,'k-',h,e2,'k--',h,e3,'k-.');
hold on;
xlabel('timestep (h)'); ylabel('numerical error');
loglog(h,error_dirk2a,'r-o',h,error_dirk3a,'b-o');
legend('slope 1','slope 2','slope 3','DIRK2','DIRK3',location='southeast');
hold off;
title(['DIRK2 and DIRK3 maximum error vs. h, y0 = ',num2str(y0_a)]);

%% part (b)
% Repeat the task with y(0) = sin(pi/4) + 10. To understand what is going
% on, plot |e(t)| for each method where e(t) is the difference between the
% numerical and exact solutions for three values of h: h = 10^-1, 10^-2,
% and 10^-3. Set the log scale in the y-axis. Do so for Tmax = 10 and Tmax
% = 1. 
y0_b = sin(pi/4) + 10;

error_dirk2b = zeros(1,length(h));
error_dirk3b = zeros(1,length(h));
for j = 1:length(h)
    [error_dirk2b(j),~] = DIRK2(gamma1,Tmax,L,h(j),y0_b);
    [error_dirk3b(j),~] = DIRK3(gamma2,Tmax,L,h(j),y0_b);
end

figure(2);
c1 = 7;
c2 = 15000000;
c3 = 300000000000;
e1 = c1*ones(1,length(h));
e2 = c2*h.^2;
e3 = c3*h.^3;
loglog(h,e1,'k-',h,e2,'k--',h,e3,'k-.');
hold on;
loglog(h,error_dirk2b,'r-o',h,error_dirk3b,'b-o');
xlabel('timestep (h)'); ylabel('numerical error');
legend('slope 0','slope 2','slope 3','DIRK2','DIRK3',location='southeast');
hold off;
title(['DIRK2 and DIRK3 maximum error vs. h, y0 = ',num2str(y0_b)]);


%% to investigate funky behavior in figure 2:
Tmax1 = 1;
h_b = [1e-1, 1e-2, 1e-3];

for k = 1:length(h_b)
    t1 = 0:h_b(k):Tmax1;
    t10 = 0:h_b(k):Tmax;

    [~,e_t10_dirk2] = DIRK2(gamma1,Tmax,L,h_b(k),y0_b);
    [~,e_t10_dirk3] = DIRK3(gamma2,Tmax,L,h_b(k),y0_b);
    [~,e_t1_dirk2] = DIRK2(gamma1,Tmax1,L,h_b(k),y0_b);
    [~,e_t1_dirk3] = DIRK3(gamma2,Tmax1,L,h_b(k),y0_b);

    % Tmax = 10 plot
    figure(3);
    hold on;
    semilogy(t10,e_t10_dirk2,t10,e_t10_dirk3);
    xlabel('time'); ylabel('|e(t)|');
    axis([0 3 0 8]); % so you can actually see what's going on near t=0
    title(['|e(t)| for Tmax = 10, y0 = ',num2str(y0_b)]);
    legend('dirk2 h=1e-1','dirk3 h=1e-1','dirk2 h=1e-2',...
        'dirk3 h=1e-2', 'dirk2 h=1e-3','dirk3 h=1e-3');

    % Tmax = 1 plot
    figure(4);
    hold on;
    semilogy(t1,e_t1_dirk2,t1,e_t1_dirk3);
    xlabel('time'); ylabel('|e(t)|');
    title(['|e(t)| for Tmax = 1, y0 = ',num2str(y0_b)]);
    legend('dirk2 h=1e-1','dirk3 h=1e-1','dirk2 h=1e-2',...
        'dirk3 h=1e-2', 'dirk2 h=1e-3','dirk3 h=1e-3');
end


% function to calculate the b(t) linear part of the implicit formula
function b = B(t,L)
    b = L * sin(t + pi/4) + cos(t + pi/4);
end

% function to calculate the exact solution to Prothero-Robinson problem
function y = exact_sol(t,L,y0)
    y = exp(-L*t)*(y0 - sin(pi/4)) + sin(t + pi/4);
end

% Function to implement DIRK2
function [e_h,e_t] = DIRK2(g,Tmax,L,h,y0)
    t = 0:h:Tmax;
    u = y0;
    A = -L;
    e_t = zeros(1,length(t));       % for part b comparison
    e_h = abs(u - exact_sol(0,L,y0));
    e_t(1) = e_h;

    for i = 1:length(t)-1
        k1 = (A * u + B(t(i)+g*h,L)) / (1 - h * g * A);
        k2 = (A * (u + h*(1-g)*k1) + B(t(i)+h,L)) / (1 - h * g * A);
        u = u + (1-g) * h * k1 + g * h * k2;
        
        y = exact_sol(t(i+1), L, y0);
        e_t(i+1) = abs(u - y);
        e_h = max([e_h, abs(u - y)]);
    end
end

% function to implement DIRK3
function [e_h,e_t] = DIRK3(g,Tmax,L,h,y0)
    t = 0:h:Tmax;
    u = y0;
    A = -L;

    e_t = zeros(1,length(t));       % for part b comparison
    e_h = abs(u - exact_sol(0,L,y0));
    e_t(1) = e_h; 

    for i = 1:length(t)-1
        k1 = (A * u + B(t(i)+g*h,L)) / (1 - h * g * A);
        k2 = (A * (u + h*(1-2*g)*k1) + B(t(i)+(1-g)*h,L)) / (1 - h * g * A);
        u = u + 0.5 * h * k1 + 0.5 * h * k2;

        y = exact_sol(t(i+1), L, y0);
        e_t(i+1) = abs(u - y);
        e_h = max([e_h, abs(u - y)]);
    end
end
