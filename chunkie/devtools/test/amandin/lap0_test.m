% Solve Laplace PDE with a sphere boundary
clearvars;
close all;

%% geometry: sphere
% only discretize the generating curve of the sphere
[chnkr] = get_sphere_geometry();
plot(chnkr,'b-x')
hold on
quiver(chnkr,'r')

% everything is in cylindrical coordinates but ignore theta part because
% its irrelevant due to axisymmtric boundary condition and solution
npts  = chnkr.npt; % the total number of points on the sphere
src = chnkr.r(:,:); % generating curve points [r; z]

%% axisymmetric boundary condition

% place a charge of strength 1 inside the sphere on the z-axis
strength = 1.0;
charge = [0.0;0.0;0.4]; % cartesion coordinates of the charge

x=src(1,:);
y=src(1,:).*0;
z=src(2,:);
pts=[x;y;z]; % cartesion coordinates of the generating curve

% set up dirichlet b.c.
r = vecnorm(pts - charge).';
g = strength ./ (4*pi*r);

% normally would have to take the fourier transform of the boundary
% condition if the boundary condition is not axisymmetric

%% Integral equation
% 0th fourier mode of double layer potential kernel
S0 = kernel('axissymlap','s',1);

% creates the matrix mapping boundary points to boundary points
opts = [];
A = chunkermat(chnkr, S0, opts);

% solve for the density
sigma = gmres(A, g, [], 1e-12, npts);

%% Evaluate at off surface point
target = [3.0;0.0;2.0]; % target in cartesion coordinates
target_cyl = [3.0;2.0]; % target in cylindrical coordinates

% evaluate at target
% chunkerkerneval evaluate the kernel at the target using adaptive 
% quadrature with our solved density sigma
opts = [];
u_sol = chunkerkerneval(chnkr, S0, sigma, target_cyl, opts);

% compute the exact solution explicitly
r = vecnorm(target - charge).';
u_true = strength./(4*pi*r);

% compute the error
error = (u_sol-u_true)/norm(u_true)


%% geometry functions
function [chnkobj] = get_sphere_geometry()
    pref = [];
    pref.k = 16; % points per chunk

    cparams = [];
    cparams.eps = 1.0e-10;
    cparams.nover = 1;
    cparams.ifclosed = false;
    cparams.ta = -pi/2;
    cparams.tb = pi/2;
    cparams.maxchunklen = 2;
    %cparams.nchmin = 8;

    narms = 0;
    amp = 0.0;

    chnkobj = chunkerfunc(@(t) starfish(t, narms, amp), cparams, pref); 
    chnkobj = sort(chnkobj);
end