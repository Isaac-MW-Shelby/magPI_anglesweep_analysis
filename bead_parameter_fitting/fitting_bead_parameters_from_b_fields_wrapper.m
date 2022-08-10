function [bead_parameters] = ...
    fitting_bead_parameters_from_b_fields_wrapper(...
    angle_sweep_location, angle_index_input, varinput)

%% function "inputs" (not changed as often as those left as real inputs)

dummy_var = varinput;

do_save = true; % save the output/figures from this fit 

% fit with a fixed z height (useful debugging tool and for large beads)
fixed_z_fitting = false; 

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
mT_to_uT = 1000; % fun to code in numbers =)

figure_tracking_number = 1; % starts tracking figure numbers

linux_computer = true; % to use slashes correctly! 

bead_z = -abs(varinput);

%% load in data 

load(angle_sweep_location, 'data');

data_file_name = data.raw_data_name;

data_file_name = data_file_name(1:end-5);

angle_index = angle_index_input;

mt_angle = data.angle_space(angle_index);

if medfilt_data
    measured_b_fields(:, :, 1) = medfilt2(squeeze(data.mean_subtracted_fields_uT(angle_index, :, :, 1))); %#ok<UNRCH>
    measured_b_fields(:, :, 2) = medfilt2(squeeze(data.mean_subtracted_fields_uT(angle_index, :, :, 2)));
    measured_b_fields(:, :, 3) = medfilt2(squeeze(data.mean_subtracted_fields_uT(angle_index, :, :, 3)));
else
    measured_b_fields = squeeze(data.mean_subtracted_fields_uT(angle_index, :, :, :));
end

applied_field_uT = squeeze(data.mean_fields_uT(angle_index, :));

%[measured_b_fields, ~] = center_Bxyz_in_uT_image(measured_b_fields);

input_sidelength = 25;
edge_cropping = 0;
[measured_b_fields, ~, xrange_of_centered, yrange_of_centered] = ...
    center_Bxyz_in_uT_image_fixed_size(measured_b_fields, ...
                                    input_sidelength, edge_cropping);

% measured_b_fields = data_to_fit;

measured_b_fields_uT = measured_b_fields;

% intermediate_data = gradient_mask_simulated_b(measured_b_fields_uT, ...
%     gradient_value, 1, step_size);

measured_b_fields = measured_b_fields / mT_to_uT;


if isfield(data, 'precisions_uT')
    measured_field_precisions_mT = ...
        data.precisions_uT(angle_index,xrange_of_centered, yrange_of_centered, :)  / mT_to_uT;
    
else
    measured_field_precisions_mT = ...
        ones(size(measured_b_fields))*0.001; %mT, 1uT is from experimental data
end



nv_spacing_input = 45; % nm

subdivisions_per_axis = round(step_size / nv_spacing_input); % number nv centers per side length
nv_spacing = step_size / subdivisions_per_axis;

pixel_size = size(measured_b_fields);
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

% figure()
% 
% subplot(2,3,1)
% imagesc((squeeze(measured_b_fields(:, :, 1)))*mT_to_uT)
% title(['centered Bx mt_angle = ' num2str(mt_angle)])
% colorbar
% colormap(linspecer)
% pbaspect([1 1 1])
% xticklabels({''})
% yticklabels({''})
% caxis([caxis_min caxis_max])
% 
% subplot(2,3,2)
% imagesc((squeeze(measured_b_fields(:, :, 2)))*mT_to_uT)
% title('centered By')
% colorbar
% colormap(linspecer)
% pbaspect([1 1 1])
% xticklabels({''})
% yticklabels({''})
% caxis([caxis_min caxis_max])
% 
% subplot(2,3,3)
% imagesc((squeeze(measured_b_fields(:, :, 3)))*mT_to_uT)
% title('centered Bz')
% colorbar
% colormap(linspecer)
% pbaspect([1 1 1])
% xticklabels({''})
% yticklabels({''})
% caxis([caxis_min caxis_max])
% 
% subplot(2,3,4)
% imagesc((squeeze(intermediate_data(:, :, 1))))
% title('Bx thresholded')
% colorbar
% colormap(linspecer)
% pbaspect([1 1 1])
% xticklabels({''})
% yticklabels({''})
% caxis([caxis_min caxis_max])
% 
% subplot(2,3,5)
% imagesc((squeeze(intermediate_data(:, :, 2))))
% title('By thresholded')
% colorbar
% colormap(linspecer)
% pbaspect([1 1 1])
% xticklabels({''})
% yticklabels({''})
% caxis([caxis_min caxis_max])
% 
% subplot(2,3,6)
% imagesc((squeeze(intermediate_data(:, :, 3))))
% title('Bz thresholded')
% colorbar
% colormap(linspecer)
% pbaspect([1 1 1])
% xticklabels({''})
% yticklabels({''})
% caxis([caxis_min caxis_max])
% 
% drawnow

%% sets up the X, Y values of the locations of the B fields

n_x = size(measured_b_fields, 1);
n_y = size(measured_b_fields, 2);

x = step_size*linspace(0,n_x-1,n_x);
y = step_size*linspace(n_y-1, 0, n_y);
[X_unrot, Y_unrot] = meshgrid(x,y);


axes_rotation = (pi/180)*(360-225); % converted to radians


X = X_unrot*cos(axes_rotation) + Y_unrot*sin(axes_rotation);
Y = X_unrot*(-sin(axes_rotation)) + Y_unrot*cos(axes_rotation);

X = X - mean(X(:));
Y = Y - mean(Y(:));

n_x_fine = subdivisions_per_axis*size(measured_b_fields, 1);
n_y_fine = subdivisions_per_axis*size(measured_b_fields, 2);

x_fine = nv_spacing*linspace(0,n_x_fine-1,n_x_fine);
y_fine = nv_spacing*linspace(n_y_fine-1, 0, n_y_fine);
[X_unrot_fine, Y_unrot_fine] = meshgrid(x_fine,y_fine);


axes_rotation = (pi/180)*(360-225); % converted to radians


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
initial_guess(6) = 50*max(measured_b_fields(:)) * abs(initial_guess(3))^3; % in mT nm^3

title_check = '020822_E1421_r1';

if length(data_file_name) == length(title_check)

    if data_file_name == title_check
        
        initial_guess(1) = -513; % nm from center of FOV
        initial_guess(2) = 185; % nm from center of FOV
        initial_guess(3) = bead_z; % nm above the surface (z=0)
        initial_guess(4) = 1.1; % radians
        initial_guess(5) = 1.9; % radians
        initial_guess(6) = 1.85*10^9; % in mT nm^3
        
    end
end

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
    measured_b_fields, measured_field_precisions_mT, rescale_vector, ...
    gradient_value, mT_to_uT, rescale_data_to_fit, ...
    pixel_side_length, subdivisions_per_axis, ...
    front_binning_matrix, back_binning_matrix, fixed_z_fitting, ...
    mask_max_dist, mask_min_dist, debug);

if sign(bead_parameter_fits(6)) < 0
    bead_parameter_fits(6) = -bead_parameter_fits(6);
    bead_parameter_fits(4) = pi - bead_parameter_fits(4);
    bead_parameter_fits(5) = mod(bead_parameter_fits(5) + pi, 2*pi);
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
    Y_fine, nvDepth, pixel_side_length, subdivisions_per_axis, gradient_value, mT_to_uT, ...
    front_binning_matrix, back_binning_matrix);

fit_bead_fields_uT = fit_bead_fields*mT_to_uT;


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

current_fig = figure(figure_tracking_number);
figure_tracking_number = figure_tracking_number + 1;
clf

current_fig.Name = [data_file_name '_measured_fit_residual_bead_fields_for_' num2str(mt_angle)];

subplot(3,3,1)
imagesc(bx_to_plot)
title(['measured Bx mt angle = ' num2str(mt_angle)])
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

%% save outputs

current_location = pwd;


if fixed_z_fitting

    modification_suffix = ['_z_val_fixed_' num2str(abs(bead_z))];
else
    modification_suffix = []; %#ok<UNRCH>
end

if linux_computer
    slash = '/'; 
else
    slash = '\';  %#ok<UNRCH>
end


if do_save
    
    fit_data.angle = data.angle_space(angle_index);
    fit_data.raw_data_file = data_file_name;
    fit_data.bead_parameter_fits = bead_parameter_fits;
    fit_data.bead_parameter_precisions = bead_parameter_fit_precisions;
    fit_data.initial_guess = initial_guess;
    fit_data.fit_units = {'nm' 'nm' 'nm' 'radian' 'radian' 'mT nm^3 / (mu0/4pi)'};
    fit_data.measured_b_fields = measured_b_fields_uT;
    fit_data.fit_bead_fields_uT = fit_bead_fields_uT;
    fit_data.applied_field_uT = applied_field_uT;
    
    
    fit_data.mask_lims = [mask_min_dist mask_max_dist];
    
    % full_save_folder_path = [prefix data_file_name slash 'bead_parameter_fits'];
    
    bead_parameter_fit_folder = 'bead_parameter_fit_outputs';
    
    if ~exist(bead_parameter_fit_folder, 'dir')
        mkdir(bead_parameter_fit_folder)
    end
    
    % save_path = [full_save_folder_path slash 'angle_' num2str(data.angle_space(angle_index)) modification_suffix ];
    
    
    save_path = [bead_parameter_fit_folder slash 'bead_parameter_fits_' data_file_name '_angle_' num2str(data.angle_space(angle_index)) modification_suffix];
    save_path = replace(save_path, ' ', '_');
    save_path = replace(save_path, '.', 'point');
    
    save(save_path, 'fit_data')
    
    bead_parameter_fit_folder_figures = 'bead_parameter_fit_figures';
    
    if ~exist(bead_parameter_fit_folder_figures, 'dir')
        
        mkdir(bead_parameter_fit_folder_figures);
        
    end
    
    cd(bead_parameter_fit_folder_figures);
    
    
    filetype_list = {'png', 'fig', 'eps'};
    
    
    for i = 1:length(filetype_list)
        
        filetype = filetype_list{i};
        
        if ~exist(filetype, 'dir')
            mkdir(filetype);
        end
        
    end
    
    cd(current_location);
    
    for i = 1:figure_tracking_number-1
        
        current_fig = figure(i);



        fig_title = current_fig.Name;
        fig_title = replace(fig_title, ' ', '_');
        
        fig_title = [data_file_name '_' fig_title]; %#ok<AGROW>
        
        for j = 1:3
            
            filetype = filetype_list{j};
        
            full_save_string = [bead_parameter_fit_folder_figures slash ...
                 filetype slash fig_title modification_suffix];
            doPageFormat(2*[5,3]);
        
            full_save_string = replace(full_save_string, '.', 'point');
        
            saveas(current_fig, full_save_string, filetype);
        end
        
    end
    
    
    
end


end