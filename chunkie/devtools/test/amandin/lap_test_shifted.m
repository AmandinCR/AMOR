% Solve Laplace PDE with a sphere boundary axisymmetric about (1,2,t)
clearvars;
close all;
format long e;

%% geometry: sphere
% only discretize the generating curve of the sphere
[chnkr] = get_sphere_geometry();
%plot(chnkr,'b-x')
%hold on
%quiver(chnkr,'r')

% everything is in cylindrical coordinates but ignore theta part because
% we only need the generating curve
npts  = chnkr.npt; % the total number of points on the sphere
src = chnkr.r(:,:); % generating curve points [r; z]

p_modes = 5; % number of positive fourier modes
n_modes = 2*p_modes + 1; % number of fourier modes (must be odd for pos/0/neg)

%% axisymmetric boundary condition

% x,y shift in the axis of symmetry
shift = [1;2];

% place a charge of strength 1 inside the sphere on the z-axis
strength = 1.0;
charge = [1.3;2.3;0.4]; % cartesion coordinates of the charge

g = zeros(n_modes,npts);
for i=1:n_modes
    % get polar/cartesian coordinates
    theta = (i-1)*2*pi/n_modes;
    
    % cartesion coordinates of the generating curve
    x=src(1,:).*cos(theta) + shift(1);
    y=src(1,:).*sin(theta) + shift(2);
    z=src(2,:);
    pts=[x;y;z]; 

    % set up dirichlet b.c.
    r = vecnorm(pts - charge);
    g(i,:) = strength ./ (4*pi*r);
end

% Reorder FFT output to match
modes = -p_modes:p_modes;
g_fft = fft(g, n_modes, 1) / n_modes; % FFT (normalized)
g_m = fftshift(g_fft, 1);  % puts negative freqs first

%% Integral equation
% need to augment our representation for existence
mu_m = zeros(n_modes,npts);
alpha_m = zeros(n_modes,1);

opts = [];
w = chnkr.wts(:);
constraint = w.';
for i = 1:n_modes
    mode = modes(i);
    m = abs(mode) + 1;
    
    % mth fourier mode of double layer potential kernel
    D = kernel('axissymlap','d',m);

    % exterior Dirichlet double-layer jump relation
    D_m = chunkermat(chnkr, D, opts) - 0.5*eye(npts);

    rhs = g_m(i,:).';
    if mode == 0
        % Complete the m=0 double layer with alpha*S[1]
        S = kernel('axissymlap','s',m);

        one_density = ones(npts,1);
        S_m = chunkermat(chnkr, S, opts);
        s_col = S_m * one_density;

        % Augmented system
        Aug = [D_m, s_col;
               constraint, 0];
        rhs_aug = [rhs; 0];

        sol = gmres(Aug, rhs_aug, [], 1e-12, npts+1);

        mu_m(i,:) = sol(1:npts).';
        alpha_m(i) = sol(end);
    else
        % Nonzero Fourier modes do not need the augmented representation
        mu_m(i,:) = gmres(D_m, rhs, [], 1e-12, npts).';
    end
end

%% Evaluate at off surface point
target = [3.0;0.4;2.0]; % target in cartesian coordinates

% target in cylindrical coordinates (r,theta,z)
target_cyl = [sqrt((target(1)-shift(1))^2 + (target(2)-shift(2))^2);
              atan2(target(2)-shift(2), target(1)-shift(1));
              target(3)];
target_new = [target_cyl(1); target_cyl(3)];

u_sol = 0;
opts = [];
for i = 1:n_modes
    mode = modes(i);
    m = abs(mode) + 1;

    D = kernel('axissymlap','d',m);
    u_m = chunkerkerneval(chnkr, D, mu_m(i,:), target_new, opts);

    if mode == 0
        S = kernel('axissymlap','s',m);
        one_density = ones(npts,1);

        u_S = chunkerkerneval(chnkr, S, one_density, target_new, opts);
        u_m = u_m + alpha_m(i)*u_S;
    end

    % Fourier recomposition
    u_sol = u_sol + u_m * exp(1i * mode * target_cyl(2));
end
u_sol = real(u_sol);

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
