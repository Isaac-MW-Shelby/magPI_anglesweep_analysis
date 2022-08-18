

%% generate B field maps from freq scan 

slash = '/';

input_signed_resonance_order = [2 -3 -1 4];

load_new_resmap = true;


if load_new_resmap

    [file, path, indx] = uigetfile;

    load([path slash file] , 'anglesweep_WFfullresonancemaps')
    
elseif ~exist('anglesweep_WFfullresonancemaps', 'var') %#ok<UNRCH>
    disp('anglesweep_WFfullresonancemaps variable not found, set load_new_resmap to true')
    return
end

[Bx, By, Bz] = generate_B_field_maps_from_freq_scan(anglesweep_WFfullresonancemaps, input_signed_resonance_order);


GHz_per_T = 28;
uT_per_GHz = (10^6)/GHz_per_T;

mean_Bx = mean(Bx(:));
mean_By = mean(By(:));
mean_Bz = mean(Bz(:));

mean_subtracted_Bx_uT = uT_per_GHz*(Bx - mean_Bx);
mean_subtracted_By_uT = uT_per_GHz*(By - mean_By);
mean_subtracted_Bz_uT = uT_per_GHz*(Bz - mean_Bz);

caxis_lim = 2; % uT


figure();

subplot(2,2,1)

imagesc(mean_subtracted_Bx_uT)
pbaspect([ 1 1 1]);
xticks([])
yticks([])
colorbar
caxis([-caxis_lim, caxis_lim])
title('Bx (\muT)')

subplot(2,2,2)

imagesc(mean_subtracted_By_uT)
pbaspect([ 1 1 1]);
xticks([])
yticks([])
colorbar
caxis([-caxis_lim, caxis_lim])
title('By (\muT)')

subplot(2,2,3)

imagesc(mean_subtracted_Bz_uT)
pbaspect([ 1 1 1]);
xticks([])
yticks([])
colorbar
caxis([-caxis_lim, caxis_lim])
title('Bz (\muT)')

disp('displaying B field components, press any key to continue to fitting')

pause()


disp('proceeding to bead fitting')

%% set up fitting for bead moment from B fields generated

do_save = true; % save the output/figures from this fit 

% fit with a fixed z height (useful debugging tool and for large beads)
fixed_z_fitting = true; 

rescale_data_to_fit = false; % rescales data and simulation -1 -> 1

debug = false; % additional figures/outputs for debugging

medfilt_data = false; % median filter data before fitting

using_turbobeads = true;

% mask out center region

mask_max_dist = 100;
mask_min_dist = 2;

bead_z = -200; % nm, for fixed z fitting 
gradient_value = 0.04; % in uT per nm, gradient at which signal lose occurs
step_size = 560; % nm, pixel
pixel_side_length = step_size;
nvDepth = 75; % nm, depth of NV layer from the surface of the diamond 
uT_per_mT = 1000; % fun to code in numbers =)

figure_tracking_number = 1; % starts tracking figure numbers


%% load in data 

measured_b_fields_uT = zeros([size(mean_subtracted_Bx_uT) 3]);

if medfilt_data
    measured_b_fields_uT(:, :, 1) = medfilt2(mean_subtracted_Bx_uT); %#ok<UNRCH>
    measured_b_fields_uT(:, :, 2) = medfilt2(mean_subtracted_By_uT);
    measured_b_fields_uT(:, :, 3) = medfilt2(mean_subtracted_Bz_uT);
else
    measured_b_fields_uT(:, :, 1) = mean_subtracted_Bx_uT;
    measured_b_fields_uT(:, :, 2) = mean_subtracted_By_uT;
    measured_b_fields_uT(:, :, 3) = mean_subtracted_Bz_uT;
end

pre_cropped_fields = measured_b_fields_uT;


mean_Bx = mean(Bx(:));
mean_By = mean(By(:));
mean_Bz = mean(Bz(:));

applied_field_uT = [mean_Bx mean_By mean_Bz];

%[measured_b_fields, ~] = center_Bxyz_in_uT_image(measured_b_fields);

input_sidelength = 25;
edge_cropping = 0;
[measured_b_fields_uT, ~, xrange_of_centered, yrange_of_centered] = ...
    center_Bxyz_in_uT_image_fixed_size(measured_b_fields_uT, ...
                                    input_sidelength, edge_cropping);
                                
 
% intermediate_data = gradient_mask_simulated_b(measured_b_fields_uT, ...
%     gradient_value, 1, step_size);

measured_b_fields_mT = measured_b_fields_uT / uT_per_mT;



measured_field_precisions_mT = ...
    ones(size(measured_b_fields_mT))*0.0006; %mT, 0.6 uT is from experimental data




nv_spacing_input = 45; % nm

subdivisions_per_axis = round(step_size / nv_spacing_input); % number nv centers per side length
nv_spacing = step_size / subdivisions_per_axis;

pixel_size = size(measured_b_fields_mT);
n_x = pixel_size(2)*subdivisions_per_axis;
n_y = pixel_size(1)*subdivisions_per_axis;


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

%% initial guess

bead_parameters = zeros(1,6);
bead_x = 0; % nm
bead_y = 0; % nm
radian_offset = 15*pi/180;
bead_theta = pi/2 - radian_offset;


% applied field from 010422_r31 at MT_theta = 0
applied =  1.0e+04 * [-0.8981    4.2156    1.1702];

% bead_phi = pi/2+pi/8;
bead_phi = atan(-applied(1) / applied(2)) + pi/2; 

bead_m = 7*10^8; % mu0/4pi*moment mag in mT nm^3
% (mu_0/4pi)*10^-15 Am^2 is equal to 10^8 mT nm^3 in these units

bead_parameters(1) = bead_x;
bead_parameters(2) = bead_y;
bead_parameters(3) = bead_z;
bead_parameters(4) = bead_theta;
bead_parameters(5) = bead_phi;
bead_parameters(6) = bead_m;

%% plot data to fit (commented out right now)



if using_turbobeads
    caxis_lim = 2;
else
    caxis_lim = 6;
end



caxis_min = -caxis_lim;
caxis_max = caxis_lim;

figure()

subplot(2,2,1)
imagesc((squeeze(measured_b_fields_uT(:, :, 1))))
title(['centered Bx'])
colorbar
colormap(linspecer)
pbaspect([1 1 1])
xticklabels({''})
yticklabels({''})
caxis([caxis_min caxis_max])

subplot(2,2,2)
imagesc((squeeze(measured_b_fields_uT(:, :, 2))))
title('centered By')
colorbar
colormap(linspecer)
pbaspect([1 1 1])
xticklabels({''})
yticklabels({''})
caxis([caxis_min caxis_max])

subplot(2,2,3)
imagesc((squeeze(measured_b_fields_uT(:, :, 3))))
title('centered Bz')
colorbar
colormap(linspecer)
pbaspect([1 1 1])
xticklabels({''})
yticklabels({''})
caxis([caxis_min caxis_max])


drawnow

%% sets up the X, Y values of the locations of the B fields

n_x = size(measured_b_fields_mT, 1);
n_y = size(measured_b_fields_mT, 2);

x = step_size*linspace(0,n_x-1,n_x);
y = step_size*linspace(n_y-1, 0, n_y);
[X_unrot, Y_unrot] = meshgrid(x,y);


axes_rotation = (pi/180)*(360-225); % converted to radians


X = X_unrot*cos(axes_rotation) + Y_unrot*sin(axes_rotation);
Y = X_unrot*(-sin(axes_rotation)) + Y_unrot*cos(axes_rotation);

X = X - mean(X(:));
Y = Y - mean(Y(:));

n_x_fine = subdivisions_per_axis*size(measured_b_fields_mT, 1);
n_y_fine = subdivisions_per_axis*size(measured_b_fields_mT, 2);

x_fine = nv_spacing*linspace(0,n_x_fine-1,n_x_fine);
y_fine = nv_spacing*linspace(n_y_fine-1, 0, n_y_fine);
[X_unrot_fine, Y_unrot_fine] = meshgrid(x_fine,y_fine);


X_fine = X_unrot_fine*cos(axes_rotation) + Y_unrot_fine*sin(axes_rotation);
Y_fine = X_unrot_fine*(-sin(axes_rotation)) + Y_unrot_fine*cos(axes_rotation);

X_fine = X_fine - mean(X_fine(:));
Y_fine = Y_fine - mean(Y_fine(:));

% figure(); imagesc(X); colorbar; pbaspect([1 1 1]); colormap(linspecer);
% title('plot of X')
% 
% figure(); imagesc(Y); colorbar; pbaspect([1 1 1]); colormap(linspecer)
% title('plot of Y')

%% set up bead parameters initial guess

initial_guess = zeros(1,6);

% right now, puts guess as center of FOV, little bit more than bead radius
% off the surface, moment scaled to match max B field generated from that
% height, moment parallel to surface of diamond

% maybe worth having different "starting" positions to go from?


initial_guess(1) = 0; % nm from center of FOV
initial_guess(2) = 0; % nm from center of FOV
initial_guess(3) = bead_z; % nm above the surface (z=0)
initial_guess(4) = pi / 2; % radians
initial_guess(5) = pi; % radians 
initial_guess(6) = 100*max(measured_b_fields_mT(:)) * abs(initial_guess(3))^3; % in mT nm^3


% initial_guess = bead_parameters;

%initial_guess = fit_data.bead_parameter_fits;

%% does the fit

x_rescale = 100;
y_rescale = 100;
z_rescale = 100;
theta_rescale = 1;
phi_rescale = 1;
m_rescale = 10^8;

rescale_vector = [x_rescale y_rescale z_rescale theta_rescale phi_rescale m_rescale];

% rescale_vector = ones(size(rescale_vector));

[bead_parameter_fits, bead_parameter_fit_precisions] = fitting_bead_parameters_from_b_fields(X_fine,Y_fine,nvDepth, initial_guess, ...
    measured_b_fields_mT, measured_field_precisions_mT, rescale_vector, ...
    gradient_value, uT_per_mT, rescale_data_to_fit, ...
    pixel_side_length, subdivisions_per_axis, ...
    front_binning_matrix, back_binning_matrix, fixed_z_fitting, ...
    mask_max_dist, mask_min_dist, debug);

% levenburg marquardt (sp?) algorithm doesnt allow bounds, this restricts
% the angles and moment magnitude accordingly after the fit

% restrict angles between 0 and 2pi
bead_parameter_fits(4) = mod(bead_parameter_fits(4), 2*pi);
bead_parameter_fits(5) = mod(bead_parameter_fits(5), 2*pi);

% restrict theta between 0 and pi
if bead_parameter_fits(4) > pi
   bead_parameter_fits(4) = 2*pi - bead_parameter_fits(4); 
   bead_parameter_fits(5) = mod(bead_parameter_fits(5)+pi, 2*pi); 
end

% restrict m to be positive
if sign(bead_parameter_fits(6)) < 0
    bead_parameter_fits(6) = -bead_parameter_fits(6);
    bead_parameter_fits(4) = pi - bead_parameter_fits(4);
    bead_parameter_fits(5) = mod(bead_parameter_fits(5)+pi, 2*pi);
end



% disp('simulation parameters')
% bead_parameters(1:5)
% disp('fits parameters')
% bead_parameter_fits(1:5)
% disp('fits precisions')
% bead_parameter_fit_precisions(1:5)
% 
% disp('simulation parameters')
% bead_parameters(6)
% disp('fits parameters')
% bead_parameter_fits(6)
% disp('fits precisions')
% bead_parameter_fit_precisions(6)


%% compare fit to real data

% [fit_bead_fields] = generateBeadFields(bead_parameter_fits, ...
%     X, Y, nv_depth);

[fit_bead_fields] = simulate_bead_B_with_gradient_mask(bead_parameter_fits, X_fine, ...
    Y_fine, nvDepth, pixel_side_length, subdivisions_per_axis, gradient_value, uT_per_mT, ...
    front_binning_matrix, back_binning_matrix);

fit_bead_fields_uT = fit_bead_fields*uT_per_mT;


x_comp_fit = squeeze(fit_bead_fields_uT(:, :, 1));
y_comp_fit = squeeze(fit_bead_fields_uT(:, :, 2));
z_comp_fit = squeeze(fit_bead_fields_uT(:, :, 3));


% figure()
% 
% subplot(1,3,1)
% imagesc(x_comp)
% title('fit Bx')
% colorbar
% colormap(linspecer)
% pbaspect([1 1 1])
% xticklabels({''})
% yticklabels({''})
% caxis([caxis_min caxis_max])
% 
% subplot(1,3,2)
% imagesc(y_comp)
% title('fit By')
% colorbar
% colormap(linspecer)
% pbaspect([1 1 1])
% xticklabels({''})
% yticklabels({''})
% caxis([caxis_min caxis_max])
% 
% subplot(1,3,3)
% imagesc(z_comp)
% title('fit Bz')
% colorbar
% colormap(linspecer)
% pbaspect([1 1 1])
% xticklabels({''})
% yticklabels({''})
% caxis([caxis_min caxis_max])

% subplot(2,3,5)
% imagesc((x_comp.^2 + y_comp.^2 + z_comp.^2).^(1/2))
% title('fit B magnitude')
% colorbar
% colormap(linspecer)
% pbaspect([1 1 1])
% xticklabels({''})
% yticklabels({''})
% 
% subplot(2,3,[3 6])
% quiver3(X_unrot,Y_unrot,zeros(size(X)), x_comp , y_comp, z_comp)
% xlim([0 20])
% ylim([0 20])
% pbaspect([1 1 1])
% xticklabels({''})
% yticklabels({''})

% view(2)

%% plot residuals 


measured_bx = squeeze(measured_b_fields_uT(:, :, 1));
measured_by = squeeze(measured_b_fields_uT(:, :, 2));
measured_bz = squeeze(measured_b_fields_uT(:, :, 3));

% bx_to_plot = medfilt2(measured_bx);
% by_to_plot = medfilt2(measured_by);
% bz_to_plot = medfilt2(measured_bz);

bx_to_plot = measured_bx;
by_to_plot = measured_by;
bz_to_plot = measured_bz;

figure()
clf

current_fig.Name = ['fit vs measured for selected angle'];

subplot(3,3,1)
imagesc(bx_to_plot)
title(['measured Bx'])
colorbar
colormap(linspecer)
pbaspect([1 1 1])
xticklabels({''})
yticklabels({''})
caxis([caxis_min caxis_max])

subplot(3,3,2)
imagesc(by_to_plot)
title('measured By')
colorbar
colormap(linspecer)
pbaspect([1 1 1])
xticklabels({''})
yticklabels({''})
caxis([caxis_min caxis_max])

subplot(3,3,3)
imagesc(bz_to_plot)
title('measured Bz')
colorbar
colormap(linspecer)
pbaspect([1 1 1])
xticklabels({''})
yticklabels({''})
caxis([caxis_min caxis_max])

subplot(3,3,7)
imagesc(bx_to_plot-x_comp_fit)
title('residual Bx')
colorbar
colormap(linspecer)
pbaspect([1 1 1])
xticklabels({''})
yticklabels({''})
caxis([caxis_min caxis_max])

subplot(3,3,8)
imagesc(by_to_plot-y_comp_fit)
title('residual By')
colorbar
colormap(linspecer)
pbaspect([1 1 1])
xticklabels({''})
yticklabels({''})
caxis([caxis_min caxis_max])

subplot(3,3,9)
imagesc(bz_to_plot-z_comp_fit)
title('residual Bz')
colorbar
colormap(linspecer)
pbaspect([1 1 1])
xticklabels({''})
yticklabels({''})
caxis([caxis_min caxis_max])


subplot(3,3,4)
imagesc(x_comp_fit)
title('fit Bx')
colorbar
colormap(linspecer)
pbaspect([1 1 1])
xticklabels({''})
yticklabels({''})
caxis([caxis_min caxis_max])

subplot(3,3,5)
imagesc(y_comp_fit)
title('fit By')
colorbar
colormap(linspecer)
pbaspect([1 1 1])
xticklabels({''})
yticklabels({''})
caxis([caxis_min caxis_max])

subplot(3,3,6)
imagesc(z_comp_fit)
title('fit Bz')
colorbar
colormap(linspecer)
pbaspect([1 1 1])
xticklabels({''})
yticklabels({''})
caxis([caxis_min caxis_max])

drawnow


%% output plot of bead moment and applied field vector
