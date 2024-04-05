function hw8_2()
% FEM numerical solution of the BVP 
% -\grad \cdot (a(x,y) \grad u ) = 0, (x,y) \in \Omega = [0,3]^2
% u(left side of square) = 0
% u(right side of square) = 1
% du/dn = 0 on top and bottom of the square

% make mesh of square with circle inside using distmesh2d
fd = @(p) drectangle(p,0,3,0,3);
theta = [0:(pi/36):2*pi]';
pfix = [1.5 + cos(theta), 1.5 + sin(theta); 3,0; 3,3; 0,3; 0,0];
[pts, tri] = distmesh2d(fd,@huniform,0.0872,[0,0;3,3],pfix);
    % used 0.0872 because that is approximate distance between points on
    % the circle

right_boundary_ind = find(round(pts(:,1),3) == 3);
left_boundary_ind = find(round(pts(:,1),3) == 0);
top_boundary_ind = find(round(pts(:,2),3) == 3);
bottom_boundary_ind = find(round(pts(:,2),3) == 0);

% pts is a N-by-2 array with coordinates of the mesh points
% tri is a Ntriag-by-3 array of indices of the triangular elements

%%
% Find boundary points with homogeneous Neumann BCs
% The number of rows in neumann is the number of boundary intervals with
% Neumann BCs

% Each row of neumann contains indices of endpoints of the corresponding
% boundary interval

neumann = [top_boundary_ind, bottom_boundary_ind];

% Find boundary points with Dirichlet BCs
% dirichlet is a column vector of indices of points with Dirichlet BCs
dirichlet = [right_boundary_ind; left_boundary_ind];


%% Start numerical solution using FEM
Npts = size(pts,1);
FreeNodes = setdiff(1:Npts,dirichlet); %mesh points with unknown values of u
A = sparse(Npts,Npts);
b = sparse(Npts,1);

%% Assembly
%% The Stiffness matrix
Ntri = size(tri,1);
centroids = zeros(Ntri,2);
a_tri = zeros(Ntri,1);      % stores the values of a for each centroid

for j = 1:Ntri % for all triangles    
    % need to compute centroid of triangles
    centroids(j,:) = [sum(pts(tri(j,:),1))/3, sum(pts(tri(j,:),2))/3];
   
    % then check if in circle by using (x-1.5)^2 + (y-1.5)^2 <=1 for (x,y)
    % being the centroid
    in_circle = (centroids(j,1)-1.5)^2 + (centroids(j,2)-1.5)^2 <= 1;
    if in_circle
        % a_tri(j) = 1.2;    % a1
        a_tri(j) = 0.8;   % for second case
    else
        a_tri(j) = 1;      % a2
    end
    
    A(tri(j,:),tri(j,:)) = A(tri(j,:),tri(j,:)) ...
      + a_tri(j) * stima3(pts(tri(j,:),:));
            % stima3 computes M = 0.5*|T_j|*G*G';
end

%% The Right-hand side, i.e., the load vector
% Volume Forces
for j = 1:Ntri
  b(tri(j,:)) = 0;  % for the case where f = 0
end

% Neumann conditions
for j = 1 : size(neumann,1)
  b(neumann(j,:)) = b(neumann(j,:)) + norm(pts(neumann(j,1),:)- ...
      pts(neumann(j,2),:)) * myg(sum(pts(neumann(j,:),:))/2)/2;
end

% Dirichlet conditions 
u = sparse(size(pts,1),1);
u(dirichlet) = myu_d(pts(dirichlet,:));
b = b - A * u;

% Computation of the solution
u(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);

% graphic representation
figure(1);
trisurf(tri,pts(:,1),pts(:,2),full(u)','facecolor','interp');
hold on;
axis ij;
colorbar;
view(2);
axis([0 3 0 3]);
set(gca,'YDir','normal');
% title('Voltage Plot, a1 = 1.2, a2 = 1'); xlabel('X'); ylabel('Y');
title('Voltage Plot, a1 = 0.8, a2 = 1'); xlabel('X'); ylabel('Y');

%% find the current density at every centroid
abs_curr_centers = zeros(Ntri,1);
for j = 1:Ntri
    x1 = pts(tri(j,1),1);
    y1 = pts(tri(j,1),2);
    x2 = pts(tri(j,2),1);
    y2 = pts(tri(j,2),2);
    x3 = pts(tri(j,3),1);
    y3 = pts(tri(j,3),2);
    u1 = u(tri(j,1));
    u2 = u(tri(j,2));
    u3 = u(tri(j,3));

    a = (u2 - u1)*(y3 - y1) - (y2 - y1)*(u3 - u1);
    b = (x2 - x1)*(u3 - u1) - (u2 - u1)*(x3 - x1);
    c = (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1);

    grad_u = [a/c b/c];

    abs_curr_centers(j) = norm(-a_tri(j) * grad_u);
end

% then recompute the absolute value of current density at triangle vertices
abs_curr_verts = zeros(Npts,1);
count_tri = zeros(Npts,1);
for j = 1:Ntri
    abs_curr_verts(tri(j,:)) = abs_curr_verts(tri(j,:))...
        + abs_curr_centers(j);
    count_tri(tri(j,:)) = count_tri(tri(j,:)) + 1;
end
abs_curr_verts = abs_curr_verts ./ count_tri;

% plot absolute value of current density
figure(2);
trisurf(tri,pts(:,1),pts(:,2),full(abs_curr_verts)','facecolor','interp');
hold on;
axis ij;
colorbar;
view(2);
axis([0 3 0 3]);
set(gca,'YDir','normal');
% title('Current Density Plot, a1 = 1.2, a2 = 1'); xlabel('X'); ylabel('Y');
title('Current Density Plot, a1 = 0.8, a2 = 1'); xlabel('X'); ylabel('Y');


end



%%
function DirichletBoundaryValue = myu_d(x)
xmin = min(x(:,1));
xmax = max(x(:,1));
midx = 0.5*(xmin + xmax);
DirichletBoundaryValue =  0.5 * (sign(x(:,1) - midx) + 1);
end

%%
function Stress = myg(x)
Stress = zeros(size(x,1),1);
end

%%
function M = stima3(vertices)
d = size(vertices,2);
G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
M = det([ones(1,d+1);vertices']) * G * G' / prod(1:d);
end
