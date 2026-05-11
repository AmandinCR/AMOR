function [sigma_m,alpha_m] = blackbox_solver(chnkr, g_m)

% chnkr - the boundary used in the solver
% g_m - the boundary condition on chnkr

npts = size(g_m,2);
n_modes = size(g_m,1);
p_modes = (n_modes-1)/2;
modes = -p_modes:p_modes;

sigma_m = zeros(n_modes,npts);
alpha_m = zeros(n_modes,1);

opts = [];
for i = 1:n_modes
    mode = modes(i);
    m = abs(mode) + 1;
    
    D = kernel('axissymlap','d',m);
    D_m = chunkermat(chnkr, D, opts) - 0.5*eye(npts);

    rhs = g_m(i,:).';
    if mode == 0
        S = kernel('axissymlap','s',m);

        one_density = ones(npts,1);
        S_m = chunkermat(chnkr, S, opts);
        s_col = S_m * one_density;
        w = chnkr.wts(:);

        Aug = [D_m, s_col;
               w.', 0];
        rhs_aug = [rhs; 0];

        sol = gmres(Aug, rhs_aug, [], 1e-12, npts+1);

        sigma_m(i,:) = sol(1:npts).';
        alpha_m(i) = sol(end);
    else
        sigma_m(i,:) = gmres(D_m, rhs, [], 1e-12, npts).';
    end
end

end