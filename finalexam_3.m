% AMSC 661 Final Exam Problem 3
% Consider the traffic model p_t + [pv(p)]_x = 0, where v(p) = 1-p^2. The
% initial density is p(0,x) = p_0(x) = 0.1 if x<0, 0.1+0.8x if 0<=x<=1, and
% 0.9 if x>1. 

%% Part (a)
% Find the equation of the characteristics and plot them starting at x0 =
% -2, -1.9, -1.8, ..., 2. 

xmin = -2;
xmax = 2;
x0 = xmin : 0.1 : xmax;
x = xmin : 0.01 : xmax;

ind_x01 = find(x0 < 0);
ind_x02 = find(x0 >= 0 & x0 <= 1);
ind_x03 = find(x0 > 1);

figure(1);
hold on;
for i = 1:length(ind_x01)
    plot(x, (x-x0(ind_x01(i)))/(1-3*(0.1)^2),'b');
end
for i = 1:length(ind_x02)
    plot(x, (x-x0(ind_x02(i)))./(1-3*(0.1+0.8*x0(ind_x02(i))).^2),'r');
end
for i = 1:length(ind_x03)
    plot(x, (x-x0(ind_x03(i)))/(1-3*(0.9)^2),'g');
end
hold off;
axis([-2 2 0 1]);
xlabel('x'); ylabel('t');
title('Characteristics for Traffic Model');

%% part (d)
% Solve the equation numerically using Godunov's method on the
% computational domain [-2,2] for time [0,10]. You must choose an
% appropriate numerical flux function. Plot solution at time t=1,2,...,10.

% just have to initialize all the stuff and do the plots :)
k = .05;        % time-step
h = 0.01;      % space-step
x = -2:h:2-h;
t = 0:k:10;
u = zeros(length(x),length(t));     % x is contained in rows, t in cols

% initial condition
ind1 = find(x < 0);
u(ind1,1) = 0.1;

ind2 = find(x >= 0 & x <= 1);
u(ind2,1) = 0.1 + 0.8*x(ind2);

ind3 = find(x > 1);
u(ind3,1) = 0.9;


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

% implement periodic boundary conditions
u(end+1,:) = u(1,:);

x = [x,x(end)+h];
tt = 0:10;
t_ind = find(t==0 | t==1 | t==2 | t==3 | t==4 | t==5 | t==6 | t== 7 ...
             | t==8 | t==9 | t==10);

figure(1);
for i = 1:length(tt)
    hold on;
    plot(x, u(:,t_ind(i)));
end
hold off;
legend('t=0','t=1','t=2','t=3','t=4','t=5',...
       't=6','t=7','t=8','t=9','t=10','Location','northwest');
title('Numerical Solution for Godunov times t=1,...,10');
xlabel('x'); ylabel('u');


 

function output = fun(u)
    output = u .* (1-u.^2);
end

function output = Godunov_flux(uL,uR)
    if uL <= uR
        x = fminbnd(@(u) fun(u), uL, uR);
        output = 0.05*fun(x);
    else
        x = fminbnd(@(u) -fun(u), uR, uL);
        output = 0.05*fun(x);
    end
end

