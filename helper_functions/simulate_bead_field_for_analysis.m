

function [Bx, By, Bz] = simulate_bead_field_for_analysis(bead_x, ...
    bead_y,bead_z, bead_theta, bead_phi, bead_m, ...
    fov_sidelength_in_pixels, max_gradient_value, num_nv_per_pixel)



bead_parameters = zeros(1,6);


bead_parameters(1) = bead_x;
bead_parameters(2) = bead_y;
bead_parameters(3) = bead_z;
bead_parameters(4) = bead_theta;
bead_parameters(5) = bead_phi;
bead_parameters(6) = bead_m;

step_size = 560; % nm
%fov_sidelength_in_pixels = 20;
nv_depth = 75; % nm
pixel_size = [fov_sidelength_in_pixels fov_sidelength_in_pixels];

% num_nv_per_pixel = 12;

subdivisions_per_axis = round(num_nv_per_pixel*step_size/500); % number nv centers per side length

conversion_to_uT = 1000; % convert from mT to uT
%max_gradient_value = 0.04; % in uT per nm


n_x = pixel_size(1)*subdivisions_per_axis;
n_y = pixel_size(2)*subdivisions_per_axis;

x = linspace(0,n_x-1,n_x);
y = linspace(n_y-1, 0, n_y);

[X_unrot, Y_unrot] = meshgrid(x,y);

axes_rotation = (pi/180)*(360-225); % converted to radians

%axes_rotation = 0;

X = X_unrot*cos(axes_rotation) + Y_unrot*sin(axes_rotation);
Y = X_unrot*(-sin(axes_rotation)) + Y_unrot*cos(axes_rotation);

X = X - mean(X(:));
Y = Y - mean(Y(:));

fovxvals = X*(step_size / subdivisions_per_axis);
fovyvals = Y*(step_size / subdivisions_per_axis);

front_binning_matrix = zeros(pixel_size(1), n_x);
back_binning_matrix = zeros(n_y, pixel_size(2));

for i = 1:size(front_binning_matrix,1)
   region_of_ones = (i-1)*subdivisions_per_axis + 1: i*subdivisions_per_axis;
   front_binning_matrix(i, region_of_ones) = ones(1, length(region_of_ones));
end

for i = 1:size(back_binning_matrix,2)
    region_of_ones = (i-1)*subdivisions_per_axis + 1: i*subdivisions_per_axis;
    back_binning_matrix( region_of_ones, i) = ones(length(region_of_ones), 1);
end

data_to_fit = simulate_bead_B_with_gradient_mask(bead_parameters, fovxvals, ...
    fovyvals, nv_depth, step_size, subdivisions_per_axis, max_gradient_value, conversion_to_uT, ...
    front_binning_matrix, back_binning_matrix);

% noise = normrnd(zeros(size(data_to_fit)), (5)*(0.001)*ones(size(data_to_fit)));

%data_to_fit = data_to_fit + noise;


Bx = squeeze(data_to_fit(:, :, 1));
By = squeeze(data_to_fit(:, :, 2));
Bz = squeeze(data_to_fit(:, :, 3));


end