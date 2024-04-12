function [A,DeltaA,h,hx,hxx,hy,hyy,rhs,exact_sol] = setup()
% Boundary function A(x,y)
A = @(x,y)2*y.*sin(pi*x);
DeltaA = @(x)-2*pi^2*x(2).*sin(pi*x(1));

% differential operator applied to B(x,y) = x(1-x)yNN(x,y,v,W,u)
h = @(x)x(1).*(1-x(1)).*x(2);
hx = @(x)(1-2*x(1)).*x(2);
hxx = @(x)-2*x(2);
hy = @(x)x(1).*(1-x(1));
hyy = @(x)0;

% right-hand side
rhs = @(x)(2-pi*pi*x(2).^2).*sin(pi*x(1));
exact_sol = @(x,y)y.^2.*sin(pi*x);
end
