% AMSC 661 Homework 5 Problem 3
% (a) Integrate the system for 10 revolutions using the implicit midpoint
% rule. Plot x and y components of your numerical solutions on the same
% xy-plane. Plot the Hamiltonian vs. time for your numerical solution. Do
% this task with time steps such that there are 100, 1000, and 10000 steps
% per period. You should generate a total of 6 figures. 

period = 9.673596609249161;
z0 = [0; 0.5; 2; 0];
nt_period = [100; 1000; 10000];

% run implicit midpoint rule for part (b) and Stoermer-Verlet for part (c)
for i = 1:3
    implicit_midpt(period,nt_period(i),z0,2*i-1);
    SV(period,nt_period(i),z0,2*i+5)
end


%% Functions
function H = ham(z)
    H = 0.5*z(1)^2 + 0.5*z(2)^2 - 1/sqrt(z(3)^2 + z(4)^2);
end

function dzdt = f(z)
% function for the RHS of the ODE
    x = z(3);
    y = z(4);
    r = sqrt(x^2+y^2);
    dzdt = [-x/r^3; -y/r^3; z(1); z(2)];
end

function Df = Df(z)
% function for the Jacobian of f(z) for part (b)
    x = z(3);
    y = z(4);
    r = sqrt(x^2+y^2);
    Df = [0, 0, (3*x^2/r^5)-(1/r^3), 3*x*y/r^5;...
          0, 0, 3*x*y/r^5, (3*y^2/r^5)-(1/r^3);...
          1, 0, 0, 0;... 
          0, 1, 0, 0];
end

function k = Newton(z,h,tol,itermax)
% do the Newton iteration to find k in implicit midpoint rule step for (b)
    % the initial k found by linearizing f is:
    k = (eye(4) + (h/2)*Df(z)) \ f(z);
    for i = 1:itermax
        Fprime = eye(4) - (h/2) * Df(z + (h/2)*k);
        k = k - Fprime \ (k - f(z + (h/2)*k));
        
        if norm(k - f(z + (h/2)*k)) < tol
            break
        end
    end
end

function [] = implicit_midpt(period,nt_pd,z0,fig_num)
% do the implicit midpoint rule to integrate the system for 10 revolutions
    tol = 1e-14;
    itermax = 20;
    h = period / nt_pd;
    Nt = nt_pd * 10;
    
    z = zeros(length(z0),Nt);
    z(:,1) = z0;
    hamiltonian = zeros(Nt,1);
    hamiltonian(1) = ham(z0);

    for i = 2:Nt
        k = Newton(z(:,i-1),h,tol,itermax);
        z(:,i) = z(:,i-1) + h*k;
        hamiltonian(i) = ham(z(:,i));
    end
    
    % plot the phase plane in x and y and the hamiltonian vs. time
    figure(fig_num);
    plot(z(3,:),z(4,:),'-');
    xlabel('x'); ylabel('y');
    title(['Midpoint: Phase Plane Plot for ',num2str(nt_pd),...
        ' steps/period']);

    figure(fig_num+1);
    t = 0 : h : (period*10)-h;
    plot(t,hamiltonian);
    xlabel('t'); ylabel('H');
    title(['Midpoint: Hamiltonian vs. time for ',num2str(nt_pd),...
        ' steps/period']);
end

% for part (c) specifically
function gradu = grad_u(xy)
    x = xy(1);
    y = xy(2);
    r = sqrt(x^2 + y^2);
    gradu = [x/r^3; y/r^3];
end

function [] = SV(period,nt_pd,z0,fig_num)
% do Stoermer-Verlet on this system
    h = period / nt_pd;
    Nt = nt_pd * 10;
    
    z = zeros(length(z0),Nt);
    z(:,1) = z0;
    hamiltonian = zeros(Nt,1);
    hamiltonian(1) = ham(z0);

    for i = 2:Nt
        uv_temp = z(1:2,i-1) - (h/2) * grad_u(z(3:4,i-1));
        
        z(3:4,i) = z(3:4,i-1) + h * uv_temp;

        z(1:2,i) = uv_temp - (h/2) * grad_u(z(3:4,i));

        hamiltonian(i) = ham(z(:,i));
    end

    figure(fig_num);
    plot(z(3,:),z(4,:),'-');
    xlabel('x'); ylabel('y');
    title(['SV: Phase Plane Plot for ',num2str(nt_pd),' steps/period']);

    figure(fig_num+1);
    t = 0 : h : (period*10)-h;
    plot(t,hamiltonian);
    xlabel('t'); ylabel('H');
    title(['SV: Hamiltonian vs. time for ',num2str(nt_pd),' steps/period']);

end

