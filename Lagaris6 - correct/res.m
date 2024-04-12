function [r,dr] = res(x,v,W,u,fun,dfun,d2fun,d3fun,d4fun)
% BVP for the Poisson equation is setup here
%
% residual functions r and their derivatives w.r.t. parameters dr 
% are evaluated in this function
%
% computer diff_operator(Psi(x)) - RHS(x)
% boundary functions
[~,DeltaA,h,hx,hxx,hy,hyy,rhs,~] = setup();

% differential operator is d^2/dx^2 + d^2/dy^2
% differential operator applied to A(x,y) = y*(2*sin(pi*x)), the bdry term
d2A = DeltaA;

% differential operator applied to B(x,y) = x(1-x)y NN(x,y,v,W,u)
h = h(x);
hx = hx(x);
hxx = hxx(x);
hy = hy(x);
hyy = hyy(x);

[f,fx,fy,fxx,fyy,df,dfx,dfy,dfxx,dfyy] = NN(x,v,W,u,fun,dfun,d2fun,d3fun,d4fun);

d2B = (hxx + hyy)*f + 2*hx*fx + 2*hy*fy + h*(fxx + fyy);

% residual r = d2A + d2B - RHS
r = d2A(x) + d2B - rhs(x);

% derivative of r w.r.t. parameters - d2A and RHS don't depend on
% parameters so only derivative of d2B with respect to parameters
dr = (hxx + hyy)*df + 2*hx*dfx + 2*hy*dfy + h*(dfxx + dfyy);
end


