% Solve Laplace PDE with a sphere boundary axisymmetric about (1,2,t)
clearvars;
close all;
format long e;

%% geometry: sphere
[chnkr] = get_sphere_geometry();
npts  = chnkr.npt;
src = chnkr.r(:,:);

p_modes = 20;
n_modes = 2*p_modes + 1;
modes = -p_modes:p_modes;

shift = [-2 0; 2 0];
n_obj = size(shift,1);

%% boundary condition

strengths = [0.4613343839
-0.5060417196
-0.05514962418
-0.756003897
-0.5729080418
0.3033469187
-0.1035991993
0.1733398391
-0.9082985758
-0.9720124025
0.4246499045
-0.07288466409]';

charges = [-1.804666322	-0.6666640445	0.4810639186;
-1.712602781	-0.5961547114	-0.3978915668;
-1.692812874	-0.7598978947	-0.04037142537;
-2.374428544	0.3859961865	0.3413980127;
-1.951853327	0.821934579	0.3408871961;
-2.556072335	0.2386290097	-0.6205199157;
2.503306288	-0.09551115377	-0.1203816102;
1.906382244	0.07086328799	0.1298788333;
1.919941238	-0.5902886207	-0.5792752814;
2.453352352	0.5179215525	0.3784860049;
1.253784272	0.2603154077	0.1308362374;
2.195745763	0.7800063349	0.2001379217]';

g_full = zeros(n_obj,n_modes,npts);
for j=1:n_obj
    g = zeros(n_modes,npts);
    for i=1:n_modes
        theta = (i-1)*2*pi/n_modes;
        x=src(1,:).*cos(theta) + shift(j,1);
        y=src(1,:).*sin(theta) + shift(j,2);
        z=src(2,:);
        pts=[x;y;z]; 
    
        for k=1:length(charges)
            r = vecnorm(pts - charges(:,k));
            g(i,:) = g(i,:) + strengths(k) ./ (4*pi*r);
        end
    end
    g_full(j,:,:) = g;
end

%%

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

plot_sol(chnkr, shift, charges, strengths, sigma1_total_m, alpha1_total_m, sigma1_total_m*0, alpha1_total_m*0)

% First correction on sphere 2: u2^(1)
d1 = blackbox_chnkeval(chnkr, chnkr, shift1, shift2, sigma1_m, alpha1_m);

bc = g2 - d1;
bc_fft = fft(bc, n_modes, 1) / n_modes;
bc_m = fftshift(bc_fft, 1);

[sigma2_m, alpha2_m] = blackbox_solver(chnkr, bc_m);

% Store accumulated sphere 2 solution
sigma2_total_m = sigma2_m;
alpha2_total_m = alpha2_m;

plot_sol(chnkr, shift, charges, strengths, sigma1_total_m, alpha1_total_m, sigma2_total_m, alpha2_total_m)

d2 = blackbox_chnkeval(chnkr, chnkr, shift2, shift1, sigma2_m, alpha2_m);

%% AMR iterations
max_iter = 1;
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
target = [3.0;-3.0;-2.0]; % target in cartesian coordinates

u1 = blackbox_pteval(chnkr, shift(1,:), target, sigma1_total_m, alpha1_total_m);
u2 = blackbox_pteval(chnkr, shift(2,:), target, sigma2_total_m, alpha2_total_m);
u_sol = u1 + u2;

u_true = 0;
for k=1:length(charges)
    r = vecnorm(target - charges(:,k));
    u_true = u_true + strengths(k) ./ (4*pi*r);
end

error = (u_sol-u_true)/norm(u_true)


%% plot
plot_sol(chnkr, shift, charges, strengths, sigma1_total_m, alpha1_total_m, sigma2_total_m, alpha2_total_m)


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
