function [f,fx,fy,fxx,fyy,df,dfx,dfy,dfxx,dfyy] = NN(x,v,W,u,fun,dfun,d2fun,d3fun,d4fun)
% f = N(x,y,par) - N(x,1,par) - N_y(x,1,par) == N1 - N2 - N3
%% derivatives of the network
z = W*x + u;
s0 = fun(z); % sigma(z)
N1 = v'*s0; % N(x,y,par)
W2 = W.*W;
s1 = dfun(z); % sigma'(z)
s2 = d2fun(z); % sigma''(z)
s3 = d3fun(z); % % sigma^(3)(z)
s4 = d4fun(z); % % sigma^(4)(z)

%% derivatives of N1
N1x = v'*(W(:,1).*s1);
N1y = v'*(W(:,2).*s1);
N1xx = v'*(W2(:,1).*s2);
N1yy = v'*(W2(:,2).*s2);

%% derivatives of N2
z1 = W*[x(1);1] + u; % y == 1
s10 = fun(z1);
s11 = dfun(z1);                 % sigma'(z1)
s12 = d2fun(z1);                % sigma''(z1)
s13 = d3fun(z1);                % sigma^(3)(z1)
s14 = d4fun(z1);                % sigma^(4)(z1)
N2 = v'*s10; % N(x,1,par)
N2x = v'*(W(:,1).*s11);
N2xx = v'*(W2(:,1).*s12);

%% derivatives of N3
N3 = v'*(W(:,2).*s11);          % N_y(x,1,par)
N3x = v'*(W(:,1).*W(:,2).*s12);
N3xx = v'*(W2(:,1).*W(:,2).*s13);

%% derivatives with respect to parameters
[nv1,nv2] = size(v); % nv2 must be 1
[nw1,nw2] = size(W);
[nu1,nu2] = size(u); % nu2 must be 1
dim = nv1 + nw1*nw2 + nu1;
%
dN1 = zeros(dim,1);
dN1x = zeros(dim,1);
dN1y = zeros(dim,1);
dN1xx = zeros(dim,1);
dN1yy = zeros(dim,1);
%
dN2 = zeros(dim,1);
dN2x = zeros(dim,1);
dN2xx = zeros(dim,1);
%
dN3 = zeros(dim,1);
dN3x = zeros(dim,1);
dN3xx = zeros(dim,1);

%% Derivatives of neural network and its frst two derivatives with respect to parameters
% dN1
dN1(1:nv1) = s0; % dN1/dv
dN1(nv1+1 : nv1+nw1*nw2) = reshape((v.*s1)*(x'),[nw1*nw2,1]); % dN1/dW
dN1(nv1+nw1*nw2+1 : end) = v.*s1; % dN1/du

% dN2
dN2(1:nv1) = s10;
dN2(nv1+1 : nv1+nw1) = (v.*s11)*x(1);
dN2(nv1+nw1+1 : nv1+nw1*nw2) = v.*s11;
dN2(nv1+nw1*nw2+1 : end) = v.*s11;

% dN3
dN3(1:nv1) = W(:,2).*s11;
dN3(nv1+1 : nv1+nw1) = (v.*W(:,2).*s12)*x(1);
dN3(nv1+nw1+1 : nv1+nw1*nw2) = v.*W(:,2).*s12 + v.*s11;
dN3(nv1+nw1*nw2+1 : end) = v.*W(:,2).*s12;

%% Derivatives for N1
%% dN1x
dN1x(1:nv1) = W(:,1).*s1; % dN1x/dv
dN1x(nv1+1 : nv1+nw1*nw2) = ...
reshape((v.*W(:,1).*s2)*(x') + (v.*s1)*[1,0],[nw1*nw2,1]);
dN1x(nv1+nw1*nw2+1 : end) = v.*W(:,1).*s2;

%% dN1y
dN1y(1:nv1) = W(:,2).*s1;
dN1y(nv1+1 : nv1+nw1*nw2) = ...
reshape((v.*W(:,2).*s2)*(x') + (v.*s1)*[0,1],[nw1*nw2,1]);
dN1y(nv1+nw1*nw2+1 : end) = v.*W(:,2).*s2;

%% dN1xx
dN1xx(1:nv1) = W2(:,1).*s2;
dN1xx(nv1+1 : nv1+nw1*nw2) = ...
reshape((v.*W2(:,1).*s3)*(x') + 2*(v.*W(:,1).*s2)*[1,0],[nw1*nw2,1]);
dN1xx(nv1+nw1*nw2+1 : end) = v.*W2(:,1).*s3;

%% dN1yy
dN1yy(1:nv1) = W2(:,2).*s2;
dN1yy(nv1+1 : nv1+nw1*nw2) = ...
reshape((v.*W2(:,2).*s3)*(x') + 2*(v.*W(:,2).*s2)*[0,1],[nw1*nw2,1]);
dN1yy(nv1+nw1*nw2+1 : end) = v.*W2(:,2).*s3;

%% Derivatives for N2
%% dN2x
dN2x(1:nv1) = W(:,1).*s11; % dN2x/dv
dN2x(nv1+1 : nv1+nw1) = (v.*W(:,1).*s12)*x(1) + v.*s11; %dN2x/dW(:,1)
dN2x(nv1+nw1 + 1 : nv1+nw1*nw2) = v.*W(:,1).*s12; %dN2x/dW(:,2)
dN2x(nv1+nw1*nw2+1 : end) = v.*W(:,1).*s12; % dN2x/du

%% dN2xx
dN2xx(1:nv1) = W2(:,1).*s12;
dN2xx(nv1+1 : nv1+nw1) = (v.*W2(:,1).*s13)*x(1) + 2*(v.*W(:,1).*s12);
dN2xx(nv1+nw1+1 : nv1+nw1*nw2) = v.*W2(:,1).*s13;
dN2xx(nv1+nw1*nw2+1 : end) = v.*W2(:,1).*s13;

%% Derivatives for N3
%% dN3x
dN3x(1:nv1) = W(:,1).*W(:,2).*s12; % dN3x/dv
dN3x(nv1+1 : nv1+nw1) = (v.*W(:,1).*W(:,2).*s13)*x(1) + v.*W(:,2).*s12; % dN3x/dW
dN3x(nv1+nw1+1 : nv1+nw1*nw2) = v.*W(:,1).*W(:,2).*s13 + v.*W(:,1).*s12; % dN3x/dW
dN3x(nv1+nw1*nw2+1 : end) = v.*W(:,1).*W(:,2).*s13; % dN3x/du

%% dN3xx
dN3xx(1:nv1) = W2(:,1).*W(:,2).*s13;
dN3xx(nv1+1 : nv1+nw1) = (v.*W2(:,1).*W(:,2).*s14)*x(1) + 2*v.*W(:,1).*W(:,2).*s13; % dN3x/dW
dN3xx(nv1+nw1+1 : nv1+nw1*nw2) = v.*W2(:,1).*W(:,2).*s14 + v.*W2(:,1).*s13; % dN3x/dW
dN3xx(nv1+nw1*nw2+1 : end) = v.*W2(:,1).*W(:,2).*s14; % dN3x/du

%% Putting derivatives of the neural network together
df = zeros(dim,1);
dfx = zeros(dim,1);
dfy = zeros(dim,1);
dfxx = zeros(dim,1);
dfyy = zeros(dim,1);
%
f = N1 - N2 - N3;
fx = N1x - N2x - N3x;
fy = N1y;
fxx = N1xx - N2xx - N3xx;
fyy = N1yy;
%
df = dN1 - dN2 - dN3;
dfx = dN1x - dN2x - dN3x;
dfxx = dN1xx - dN2xx - dN3xx;
dfy = dN1y;
dfyy = dN1yy;
end