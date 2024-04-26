% AMSC 661 Homework 1e
% Set phi = max(1 - |x|,0), psi = 0, a = sqrt(2), h = 0.05, and a
% reasonable time step k. Pick the numerical domain -6 to 6 and periodic
% boundary conditions. Solve the equations for xi and eta numerically using
% Lax-Friedrichs, appropriate Upwind (left for one and right for the other
% depending on the signs of lambda), Lax-Wendroff, and appropriate
% Beam-Warming methods. Return to the variable w and then to the original
% variable u. Plot the numerical solutions u obtained using each of the
% methods at times t = 1/2a, 1/a, 2/a, and 4/a as well as the exact
% solution at these times (d'Alembert's solution). 

% solving xi_t = a*xi_x and eta_t = -a*eta_x 

method = 4;

a = sqrt(2);
h = 0.05;
k = h/a;       % because makes |ak/h|<=1

L1 = a;
L2 = -a;

t = 0:k:(4/a);
x = -6:h:(6-h);

% for each timestep, modify column of xi and eta
xi = zeros(length(x),length(t));
eta = zeros(length(x),length(t));

% initial conditions for xi and eta: 1/2 phi'(x) (weak derivative)
% at 0, we have the average of -1/2 and 1/2 (the jump discontinuities),
% which is just 0, so we don't include 0 in the modifications
% at -1 and 1, we have 0.25 and -0.25, respectively
neg_ind = find(-1<=x & x<0);
pos_ind = find(0<x & x<=1);
xi(neg_ind,1) = [0.25; 0.5*ones(length(neg_ind)-1,1)];
eta(neg_ind,1) = [0.25; 0.5*ones(length(neg_ind)-1,1)];
xi(pos_ind,1) = [-0.5*ones(length(pos_ind)-1,1); -0.25];
eta(pos_ind,1) = [-0.5*ones(length(pos_ind)-1,1); -0.25];

% making all the matrices we need
I = eye(length(x));

e = ones(length(x),1);
A1 = spdiags([-e,0*e,e],-1:1,length(x),length(x));
A1(1,end) = -1; A1(end,1) = 1;

A2 = spdiags([e,0*e,e],-1:1,length(x),length(x));
A2(1,end) = 1; A2(end,1) = 1;

A3 = spdiags([-1*e,e],-1:0,length(x),length(x));
A3(1,end) = -1;

A4 = spdiags([-1*e,e],0:1,length(x),length(x));
A4(end,1) = 1;

% method 1 - Lax-Friedrichs
if method == 1
    for i = 2:length(t)
        xi(:,i) = (0.5 * A2 - (0.5*L1*k/h) * A1) * xi(:,i-1);
        eta(:,i) = (0.5 * A2 - (0.5*L2*k/h) * A1) * eta(:,i-1);
    end
end

% method 2 - Upwind. Use left for xi (a>0) and right for eta (a<0)
if method == 2
    for i = 2:length(t)
        xi(:,i) = (I - (L1*k / h) * A3) * xi(:,i-1);
        eta(:,i) = (I - (L2*k / h) * A4) * eta(:,i-1);
    end
end

% method 3: Lax-Wendroff
if method == 3
    for i = 2:length(t)
        xi(:,i) = (I - (0.5*L1*k / h) * A1 + 0.5*(L1^2 * k^2 / h^2) ...
                   * (A2 - 2*I)) * xi(:,i-1);
        eta(:,i) = (I - (0.5*L2*k / h) * A1 + 0.5*(L2^2 * k^2 / h^2) ...
                    * (A2 - 2*I)) * eta(:,i-1);
    end
end

% method 4: Beam-Warning (left for xi and right for eta)
if method == 4
    A5 = spdiags([e,-4*e,3*e,e,-4*e],[-2,-1,0,length(x)-2,length(x)-1],...
                  length(x),length(x));
    A6 = spdiags([e,-2*e,e,e,-2*e],[-2,-1,0,length(x)-2,length(x)-1],...
                  length(x),length(x));
    
    A7 = spdiags([4*e,-e,-3*e,4*e,-e],[-length(x)+2,-length(x)+1,0,1,2],...
                  length(x),length(x));
    A8 = spdiags([-2*e,e,e,-2*e,e],[-length(x)+2,-length(x)+1,0,1,2],...
                  length(x),length(x));
    for i = 2:length(t)
        xi(:,i) = (I - 0.5*(L1*k/h)*A5 + (L1^2 * k^2 / h^2)*A6)...
                   * xi(:,i-1);
        eta(:,i) = (I - 0.5*(L2*k/h)*A7 + (L2^2 * k^2 / h^2)*A8)...
                    * eta(:,i-1);
    end
end

u_t = L1*xi + L2*eta;
u_x = xi + eta;
% use periodic boundary conditions to add the solution at x=6 to the matrix
x = [x,6];
u_x(end+1,:) = u_x(1,:);

% I choose to do numerical integration over x, at each timestep (that way I
% don't have to worry about setting the initial condition)
% use cumtrapz on xi, which will integrate each column (timestep) of xi 
u = cumtrapz(x,u_x);

% now to plot everything
sol_exact = zeros(length(x),1);

t_ind = find(t==0.5/a | t==1/a | t==2/a | t==4/a);

figure(1);
hold on;
for j = 1:length(t_ind)
    plot(x,u(:,t_ind(j)),'-','LineWidth',3);

    for k = 1:length(x)
        sol_exact(k) = exactsol(a,x(k),t(t_ind(j)));
    end

    plot(x,sol_exact,'-.','LineWidth',1);
end
xlabel('x'); ylabel('u'); 
axis([-6 6 -0.1 1]);
title(['Wave Eq. Solutions for method ',num2str(method)]);
legend('num sol t=1/2a','num sol t=1/a','num sol t=2/a','num sol t=4/a',...
    'exact t=1/2a','exact t=1/a','exact t=2/a','exact t=4/a');

function sol = exactsol(a,x,t)
    sol = 0.5 * (max(1-abs(x+a*t),0) + max(1-abs(x-a*t),0));
end





