function [b_fields] = generateBeadFields(bead_params, ...
    fov_xvals, fov_yvals, nvDepth)



x = bead_params(1);
y = bead_params(2);
z = bead_params(3);
theta = bead_params(4);
phi = bead_params(5);
moment = bead_params(6);

r_vecs = zeros([size(fov_xvals) 3]);

r_vecs(:, :, 1) = fov_xvals - x;
r_vecs(:, :, 2) = fov_yvals - y;
r_vecs(:, :, 3) = ones(size(r_vecs(:, :, 3))).*(nvDepth - z);

rmags = zeros(size(fov_xvals));

for i = 1:size(fov_xvals,1)
    for j = 1:size(fov_xvals, 2)
       rmags(i, j) = norm(squeeze(r_vecs(i,j,:))); 
    end
end

rhat = zeros(size(r_vecs));
rhat(:, :, 1) = squeeze(r_vecs(:, :, 1)) ./ rmags;
rhat(:, :, 2) = squeeze(r_vecs(:, :, 2)) ./ rmags;
rhat(:, :, 3) = squeeze(r_vecs(:, :, 3)) ./ rmags;

moment_hat = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];

b_fields = zeros([size(fov_xvals) 3]);
                       
moment_hat_dot_r_hat = squeeze(moment_hat(1).*rhat(:, :, 1) + ...
                           moment_hat(2).*rhat(:, :, 2) + ...
                           moment_hat(3).*rhat(:, :, 3));

b_fields(:, :, 1) = (moment ./ rmags.^3).*(3*moment_hat_dot_r_hat.*squeeze(rhat(:, :, 1)) - moment_hat(1));
b_fields(:, :, 2) = (moment ./ rmags.^3).*(3*moment_hat_dot_r_hat.*squeeze(rhat(:, :, 2)) - moment_hat(2));
b_fields(:, :, 3) = (moment ./ rmags.^3).*(3*moment_hat_dot_r_hat.*squeeze(rhat(:, :, 3)) - moment_hat(3));


