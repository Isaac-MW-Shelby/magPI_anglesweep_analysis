function [masked_b_field] = gradient_mask_simulated_b(simulated_b, ...
    gradient_threshold, conversion_to_uT, step_size_nm)

bx = squeeze(simulated_b(:, :, 1));
by = squeeze(simulated_b(:, :, 2));
bz = squeeze(simulated_b(:, :, 3));

normalization = 1/(sqrt(3));

o1 = normalization * [1 1 1];
o2 = normalization * [-1 -1 1];
o3 = normalization * [-1 1 -1];
o4 = normalization * [1 -1 -1];
orientations = [o1; o2 ;o3; o4];

Bprojections = zeros([size(bx) 4]);

for j = 1:4
    nhat = squeeze(orientations(j, :));
    Bprojections(:, :, j) = bx*nhat(1) + by*nhat(2) + bz*nhat(3);
end


gradient_masks = zeros([size(bx) 4]);

for j = 1:4
    [delBx,delBy] = gradient(squeeze(Bprojections(:, :, j)));
    gradient_magnitude = (delBx.^2 + delBy.^2).^(1/2);
    magnitude_grad_uTpernm = gradient_magnitude*conversion_to_uT / (step_size_nm);
    magnitude_gradmask = magnitude_grad_uTpernm < gradient_threshold;
    gradient_masks(:, :, j) = magnitude_gradmask;
end

combined_mask = ones(size(gradient_masks(:, :, 1)));

for j = 1:4
    combined_mask = combined_mask .* gradient_masks(:, :, j);
end

Bx = combined_mask.*bx;
By = combined_mask.*by;
Bz = combined_mask.*bz;

masked_b_field = zeros(size(simulated_b));

masked_b_field(:, :, 1) = Bx;
masked_b_field(:, :, 2) = By;
masked_b_field(:, :, 3) = Bz;
   


end