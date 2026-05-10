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
D0 = kernel('axissymlap','d',1);
S0 = kernel('axissymlap','s',1);

opts = [];
Dmat = chunkermat(chnkr, D0, opts);
A = Dmat - 0.5*eye(npts);

% Single-layer potential of constant density 1 on the boundary
one_density = ones(npts,1);
s_col = chunkermat(chnkr, S0, opts) * one_density;

% Add one constraint to fix the nullspace of the double layer density.
% A simple choice is mean(mu)=0.
w = chnkr.wts(:);
constraint = w.';

% Augmented system:
Aug = [A, s_col;
       constraint, 0];

rhs = [g; 0];
sol = gmres(Aug, rhs, [], 1e-12, npts+1);

mu = sol(1:npts);
alpha = sol(end);

%% Evaluate at off surface point
target = [3.0;0.0;2.0];
target_cyl = [3.0;2.0];

opts = [];

u_D = chunkerkerneval(chnkr, D0, mu, target_cyl, opts);
u_S = chunkerkerneval(chnkr, S0, one_density, target_cyl, opts);

u_sol = u_D + alpha*u_S

% compute the exact solution explicitly
r = vecnorm(target - charge).';
u_true = strength./(4*pi*r)

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