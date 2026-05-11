function [u] = blackbox_eval(chnkr, shift, target, sigma_m, alpha_m)

% chnkr - the boundary we have the solved density on
% shift - axis of symmetry for chnkr
% target - x,y,z coord of target
% sigma_m - the solved density
% alpha_m - the solved extra parameter

npts = size(sigma_m,2);
n_modes = size(sigma_m,1);
p_modes = (n_modes-1)/2;
modes = -p_modes:p_modes;

% target in cylindrical coordinates (r,theta,z)
target_cyl = [sqrt((target(1)-shift(1))^2 + (target(2)-shift(2))^2);
              atan2(target(2)-shift(2), target(1)-shift(1));
              target(3)];
target_new = [target_cyl(1); target_cyl(3)];

opts = [];
u_sol = 0;
for i = 1:n_modes
    mode = modes(i);
    m = abs(mode) + 1;

    D = kernel('axissymlap','d',m);
    u_m = chunkerkerneval(chnkr, D, sigma_m(i,:), target_new, opts);

    if mode == 0
        S = kernel('axissymlap','s',m);
        one_density = ones(npts,1);

        u_S = chunkerkerneval(chnkr, S, one_density, target_new, opts);
        u_m = u_m + alpha_m(i)*u_S;
    end

    u_sol = u_sol + u_m * exp(1i * mode * target_cyl(2));
end
u = real(u_sol);

end