% AMSC 661 Homework 3 Problem 4

% Consider stiff Robertson's problem from chemical kinetics: x' = -ax +
% byz, y' = ax - byz - cy^2, z' = cy^2, x(0) = 1, y(0) = 0, z(0) = 0. 

% (a) implement DIRKo3 from Problem 2, with the same Butcher's array and
% gamma = 0.5 + sqrt(3)/6. 

% (b) Implement the two-step BDF method. Use DIRK2 for the first time step.

% (c) Compute and plot the numerical solutions on t in [0,100] for all
% three methods at h = 10^-3, 10^-2, and 10^-1. Measure CPUtime in all
% cases. Also, plot the CPUtime vs. h for all three methods. Use log scale
% for the CPU time and for h in Figure 4. Compare the performance of these
% three methods and write a summary of your observations. 

h_vec = [1.0e-1, 1.0e-2, 1.0e-3];
Tmax = 100.0;
Nsteps = ceil(Tmax * h_vec.^(-1));
tol = 1.0e-14;
itermax = 20;

y0 = [1.0; 0.0; 0.0];

% For each value of h, find the solution for DIRK2, DIRKo3, and BDF2 and
% track their CPU times. Then plot for each value of h. 
for i = 1:length(h_vec)
    t = h_vec(i)*(1:(Nsteps(i)+1))';

    DIRK2_sol = zeros(Nsteps(i)+1,3);
    DIRK2_sol(1,:) = y0;

    DIRK3_sol = zeros(Nsteps(i)+1,3);
    DIRK3_sol(1,:) = y0;

    tic     % measure CPU time for DIRK2
    for j = 1:Nsteps(i)
        DIRK2_sol(j+1,:) = DIRK2step(DIRK2_sol(j,:)',h_vec(i),tol,itermax);
    end
    DIRK2_time = toc;

    tic     % measure CPU time for DIRK3
    for j = 1:Nsteps(i)
        DIRK3_sol(j+1,:) = DIRKo3_step(DIRK3_sol(j,:)',h_vec(i),tol,itermax);
    end
    DIRK3_time = toc;

    BDF2_sol = zeros(Nsteps(i)+1,3);
    BDF2_sol(1,:) = y0;
    % get second timestep using DIRK2 since BDF is 2-step
    BDF2_sol(2,:) = DIRK2step(y0,h_vec(i),tol,itermax);

    tic     % measure CPU time for BDF2
    for j = 2:Nsteps(i)
        BDF2_sol(j+1,:) = BDF2_step(BDF2_sol(j,:)',BDF2_sol(j-1,:)',...
                                    h_vec(i),tol,itermax);
    end
    BDF2_time = toc;


    figure(i);
    % for each h value, plot x, y, and z vs. t, all solvers on same plots
    subplot(3,1,1);
    plot(t,DIRK2_sol(:,1),'r-');
    hold on;
    plot(t,DIRK3_sol(:,1),'b--');
    plot(t,BDF2_sol(:,1),'g-.');
    hold off;
    xlabel('t'); ylabel('x');
    legend('DIRK2','DIRKo3','BDF2');
    title(['Concentrations of X, Y, and Z, h=',num2str(h_vec(i))]);

    subplot(3,1,2);
    plot(t,DIRK2_sol(:,2),'r-');
    hold on;
    plot(t,DIRK3_sol(:,2),'b--');
    plot(t,BDF2_sol(:,2),'g-.');
    hold off;
    xlabel('t'); ylabel('y');
    legend('DIRK2','DIRKo3','BDF2');

    subplot(3,1,3);
    plot(t,DIRK2_sol(:,3),'r-');
    hold on;
    plot(t,DIRK3_sol(:,3),'b--');
    plot(t,BDF2_sol(:,3),'g-.');
    xlabel('t'); ylabel('z');   
    hold off;
    legend('DIRK2','DIRKo3','BDF2');

    figure(length(h_vec)+1);
    loglog(h_vec(i),DIRK2_time,'ro');
    hold on;
    loglog(h_vec(i),DIRK3_time,'bo');
    loglog(h_vec(i),BDF2_time,'go');
    xlabel('timestep (h)'); ylabel('CPU time');
    title('CPU time vs. timestep for all methods');
    legend('DIRK2','DIRKo3','BDF2');

end
hold off;

%% BDF method
% I derived the formula F = u_n+1 - (4/3)u_n + (1/3)u_n-1 = (2/3)h*f(u_n+1)
% from the 2-step BDF derivation since h is now constant, and we will use
% Newton's method within the step solver to solve for u_n+1 implicitly. I
% also derived that DF(u_n+1) = I - (2/3)*h*Df(u_n+1), used in line 150.

function F = F_BDF(u_next,u,u_prev,h)
    F = u_next - (4/3)*u + (1/3)*u_prev - (2/3)*h*func(u_next);
end

function ynew = BDF2_step(y,y_prev,h,tol,itermax)
    ynew = y;

    for j = 1:itermax
        DF = eye(3) - (2/3)*h*Jac(ynew);
        ynew = ynew - DF \ F_BDF(ynew,y,y_prev,h);
        
        if norm(F_BDF(ynew,y,y_prev,h)) < tol
            break
        end
    end
end

%% DIRKo3
% note: solving for k1 is the same as for DIRK2, so we can use
% NewtonIterDIRK2 function for that. Need to write a new function for the
% DIRK3step which incorporates the changes for k2 and y_n+1

function ynew = DIRKo3_step(y,h,tol,itermax)
    gamma = 0.5 + sqrt(3)/6;
    k1 = func(y);
    
    % iterate until tolerance met to get k1
    for j = 1:itermax
        k1 = NewtonIterDIRK2(y,h,k1,gamma);
        
        if norm(k1 - func(y + h*gamma*k1)) < tol
            break
        end
    end
    
    % find k2 now, start it at k1
    k2 = k1;
    aux1 = y + h*(1 - 2*gamma)*k1;
    
    % iterate using Newton's method with aux1
    for j = 1:itermax
        k2 = NewtonIterDIRK2(aux1,h,k2,gamma);
        aux2 = aux1 + h*gamma*k2;
        
        if norm(k2 - func(aux2)) < tol
            break
        end
    end
    
    % Calculate next y using b' = [0.5, 0.5]
    ynew = y + 0.5*h*k1 + 0.5*h*k2;
end

%% DIRK2 (from Dr. Cameron)
function knew = NewtonIterDIRK2(y,h,k,gamma)
    aux = y + h*gamma*k;
    F = k - func(aux);
    DF = eye(3) - h*gamma*Jac(aux);
    knew = k - DF\F;
end

function ynew = DIRK2step(y,h,tol,itermax)
    gamma = 1.0 - 1.0/sqrt(2);
    k1 = func(y);
    for j = 1 : itermax
        k1 = NewtonIterDIRK2(y,h,k1,gamma);
        if norm(k1 - func(y + h*gamma*k1)) < tol
            break
        end
    end
    k2 = k1;
    y = y + h*(1-gamma)*k1;
    for j =1 : itermax
        k2 = NewtonIterDIRK2(y,h,k2,gamma);
        aux = y + h*gamma*k2;
        if norm(k2 - func(aux)) < tol
            break
        end
    end
    ynew = aux;
end

%% RHS and Jacobian functions (from Dr. Cameron)
% the right-hand side
function dy = func(y) 
    a = 0.04;
    b = 1.0e4;
    c = 3.0e7;
    dy = zeros(3,1);
    byz = b*y(2)*y(3);
    cy2 = c*y(2)*y(2);
    ax = a*y(1);
    dy(1) = -ax + byz;
    dy(2) = ax - byz - cy2;
    dy(3) = cy2;
end

% the Jacobian matrix for the right-hand side
function J = Jac(y)
    a = 0.04;
    b = 1.0e4;
    c = 3.0e7;
    by = b*y(2);
    bz = b*y(3);
    c2y = 2*c*y(2);
    J = zeros(3);
    J(1,1) = -a;
    J(1,2) = bz;
    J(1,3) = by;
    J(2,1) = a;
    J(2,2) = -bz-c2y;
    J(2,3) = -by;
    J(3,2) = c2y;
end

