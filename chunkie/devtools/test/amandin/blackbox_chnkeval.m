function [u] = blackbox_eval(chnkr, chnkr_targ, shift, shift_targ, sigma_m, alpha_m)

% chnkr - the boundary we have the solved density on
% chnkr_targ - the boundary we want the solution at
% shift - axis of symmetry for chnkr
% shift_targ - axis of symmetry for target boundary
% sigma_m - the solved density
% alpha_m - the solved extra parameter


npts = size(sigma_m,2);
n_modes = size(sigma_m,1);
p_modes = (n_modes-1)/2;
modes = -p_modes:p_modes;

src = chnkr_targ.r(:,:);
u = zeros(n_modes,npts);
for i=1:n_modes
    theta = (i-1)*2*pi/n_modes;
    x=src(1,:).*cos(theta) + shift_targ(1);
    y=src(1,:).*sin(theta) + shift_targ(2);
    z=src(2,:);
    targets = [x;y;z];
    targets_cyl = [sqrt((targets(1,:)-shift(1)).^2 + (targets(2,:)-shift(2)).^2);
                   atan2(targets(2,:)-shift(2), targets(1,:)-shift(1));
                   targets(3,:)];
    targets_new = [targets_cyl(1,:); targets_cyl(3,:)];

    u_sol = zeros(npts,1);
    opts = [];
    for k = 1:n_modes
        mode = modes(k);
        m = abs(mode) + 1;
    
        D = kernel('axissymlap','d',m);
        u_m = chunkerkerneval(chnkr, D, sigma_m(k,:), targets_new, opts);
    
        if mode == 0
            S = kernel('axissymlap','s',m);
            one_density = ones(npts,1);
            u_S = chunkerkerneval(chnkr, S, one_density, targets_new, opts);
            u_m = u_m + alpha_m(k)*u_S;
        end
        
        u_sol = u_sol + u_m .* exp(1i * mode * targets_cyl(2,:).');
    end
    u(i,:) = real(u_sol).';
end










% I want to return a n_modes x npts matrix

end