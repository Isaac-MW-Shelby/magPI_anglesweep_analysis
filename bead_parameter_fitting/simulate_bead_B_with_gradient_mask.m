
function gradient_masked_B = simulate_bead_B_with_gradient_mask(bead_parameters, fovxvals, ...
    fovyvals, nvDepth, pixel_side_length, subdivisions_per_axis, gradient_value, conversion_to_uT, ...
    front_binning_matrix, back_binning_matrix)

nv_spacing = pixel_side_length / subdivisions_per_axis;

simulated_bead_fields = generateBeadFields(bead_parameters, ...
    fovxvals, fovyvals, nvDepth);

intermediate_data = gradient_mask_simulated_b(simulated_bead_fields, ...
    gradient_value, conversion_to_uT, nv_spacing);

intermediate_data = imgaussfilt(intermediate_data, subdivisions_per_axis);

rescale_factor = (size(front_binning_matrix, 2) / size(front_binning_matrix, 1))^2; 

binned = zeros(size(front_binning_matrix, 1), size(back_binning_matrix, 2), 3);

for i = 1:3
   binned(: , :, i) =  front_binning_matrix*intermediate_data(:, :, i)*back_binning_matrix;
end

gradient_masked_B = binned / rescale_factor;

end