function plot_sol(chnkr, shift, charges, strengths, sigma1_m, alpha1_m, sigma2_m, alpha2_m)


% Grid bounds
xmin = -4; xmax = 4;
ymin = -4; ymax = 4;
Nx = 200;
Ny = 200;
xvec = linspace(xmin, xmax, Nx);
yvec = linspace(ymin, ymax, Ny);
[X,Y] = meshgrid(xvec, yvec);
Z = zeros(size(X));
targets = [X(:).'; Y(:).'; Z(:).'];

% Remove points inside the spheres if needed
a = 1.0;   % change this to your sphere radius

inside1 = sqrt((targets(1,:) - shift(1,1)).^2 + ...
               (targets(2,:) - shift(1,2)).^2 + ...
               targets(3,:).^2) < a;

inside2 = sqrt((targets(1,:) - shift(2,1)).^2 + ...
               (targets(2,:) - shift(2,2)).^2 + ...
               targets(3,:).^2) < a;

valid = ~(inside1 | inside2);

% Evaluate AMR solution
%{
u_sol_vec = nan(1, size(targets,2));
u1 = blackbox_ptseval(chnkr, shift(1,:), targets(:,valid), ...
                     sigma1_m, alpha1_m);
u2 = blackbox_ptseval(chnkr, shift(2,:), targets(:,valid), ...
                     sigma2_m, alpha2_m);
u_sol_vec(valid) = u1 + u2;
%}


% Evaluate exact solution
u_true_vec = nan(1, size(targets,2));
u_true_tmp = zeros(1, sum(valid));
for k = 1:length(strengths)
    r = vecnorm(targets(:,valid) - charges(:,k));
    u_true_tmp = u_true_tmp + strengths(k) ./ (4*pi*r);
end
u_true_vec(valid) = u_true_tmp;

%{
% Error
abs_err_vec = abs(u_sol_vec - u_true_vec);
rel_err_vec = abs(u_sol_vec - u_true_vec) ./ max(abs(u_true_vec), 1e-14);
abs_err = reshape(abs_err_vec, Ny, Nx);
rel_err = reshape(rel_err_vec, Ny, Nx);
%}


%{
% Plot absolute error
figure;
imagesc(xvec, yvec, abs_err);
axis equal tight;
set(gca,'YDir','normal');
colorbar;
title('absolute error on z = 0 slice');
xlabel('x');
ylabel('y');
%}

%{
u_sol_plot = reshape(u_sol_vec, Ny, Nx);
figure;
imagesc(xvec, yvec, u_sol_plot);
axis equal tight;
set(gca,'YDir','normal');
colorbar;
title('Approximate solution on z = 0 slice');
xlabel('x');
ylabel('y');
%}


u_true_plot = reshape(u_true_vec, Ny, Nx);
figure;
imagesc(xvec, yvec, u_true_plot);
axis equal tight;
set(gca,'YDir','normal');
colorbar;
title('True solution on z = 0 slice');
xlabel('x');
ylabel('y');


end