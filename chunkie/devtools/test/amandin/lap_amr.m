% Solve Laplace PDE with a sphere boundary axisymmetric about (1,2,t)
clearvars;
close all;
format long e;

%% geometry: sphere
[chnkr] = get_sphere_geometry();
npts  = chnkr.npt;
src = chnkr.r(:,:);

p_modes = 10;
n_modes = 2*p_modes + 1;
modes = -p_modes:p_modes;

shift = [1 2; -1 1];
n_obj = size(shift,1);

%% boundary condition

strength1 = 1.0;
charge1 = [1.5;2.5;0.4];

strength2 = 1.0;
charge2 = [-1.5;1.5;-0.4];

g_full = zeros(n_obj,n_modes,npts);
for j=1:n_obj
    g = zeros(n_modes,npts);
    for i=1:n_modes
        theta = (i-1)*2*pi/n_modes;
        x=src(1,:).*cos(theta) + shift(j,1);
        y=src(1,:).*sin(theta) + shift(j,2);
        z=src(2,:);
        pts=[x;y;z]; 
    
        r1 = vecnorm(pts - charge1);
        r2 = vecnorm(pts - charge2);
        g(i,:) = strength1 ./ (4*pi*r1) + strength2 ./ (4*pi*r2);
    end
    g_full(j,:,:) = g;
end


% boundary data on sphere1 and sphere2
g1 = reshape(g_full(1,:,:), n_modes, npts);
g2 = reshape(g_full(2,:,:), n_modes, npts);
shift1 = shift(1,:);
shift2 = shift(2,:);


% First correction on sphere 1: u1^(1)
bc = g1;
bc_fft = fft(bc, n_modes, 1) / n_modes;
bc_m = fftshift(bc_fft, 1);

[sigma1_m, alpha1_m] = blackbox_solver(chnkr, bc_m);

% Store accumulated sphere 1 solution
sigma1_total_m = sigma1_m;
alpha1_total_m = alpha1_m;

% First correction on sphere 2: u2^(1)
d1 = blackbox_chnkeval(chnkr, chnkr, shift1, shift2, sigma1_m, alpha1_m);

bc = g2 - d1;
bc_fft = fft(bc, n_modes, 1) / n_modes;
bc_m = fftshift(bc_fft, 1);

[sigma2_m, alpha2_m] = blackbox_solver(chnkr, bc_m);

% Store accumulated sphere 2 solution
sigma2_total_m = sigma2_m;
alpha2_total_m = alpha2_m;

d2 = blackbox_chnkeval(chnkr, chnkr, shift2, shift1, sigma2_m, alpha2_m);

%% AMR iterations
max_iter = 3;
for i=1:max_iter
    % sphere 1 correction
    bc = -d2;
    bc_fft = fft(bc, n_modes, 1) / n_modes;
    bc_m = fftshift(bc_fft, 1);
    
    [sigma1_m, alpha1_m] = blackbox_solver(chnkr, bc_m);

    sigma1_total_m = sigma1_total_m + sigma1_m;
    alpha1_total_m = alpha1_total_m + alpha1_m;

    d1 = blackbox_chnkeval(chnkr, chnkr, shift1, shift2, sigma1_m, alpha1_m);
    
    % sphere 2 correction
    bc = -d1;
    bc_fft = fft(bc, n_modes, 1) / n_modes;
    bc_m = fftshift(bc_fft, 1);
    
    [sigma2_m, alpha2_m] = blackbox_solver(chnkr, bc_m);

    sigma2_total_m = sigma2_total_m + sigma2_m;
    alpha2_total_m = alpha2_total_m + alpha2_m;

    d2 = blackbox_chnkeval(chnkr, chnkr, shift2, shift1, sigma2_m, alpha2_m);
end


%% error
target = [3.0;3.0;-2.0]; % target in cartesian coordinates

u1 = blackbox_pteval(chnkr, shift(1,:), target, sigma1_total_m, alpha1_total_m);
u2 = blackbox_pteval(chnkr, shift(2,:), target, sigma2_total_m, alpha2_total_m);
u_sol = u1 + u2;

r1 = vecnorm(target - charge1).';
r2 = vecnorm(target - charge2).';
u_true = strength1 ./ (4*pi*r1) + strength2 ./ (4*pi*r2);

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
