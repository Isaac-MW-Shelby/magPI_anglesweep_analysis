

%% starts by combining relevant information into a single file

current_folder_full_path = pwd;




location_of_slashes = find(current_folder_full_path == '/');

location_of_last_slash = location_of_slashes(end);

current_folder = current_folder_full_path(location_of_last_slash+1:end);

full_sweepfit_filename = ['fits_for_anglesweep_' current_folder '.mat'];

if ~exist(full_sweepfit_filename, 'file')
    
    cd('bead_parameter_fit_outputs');
    
    list_of_files = dir('*.mat');
    
    num_files = length(list_of_files);
    
    angles = zeros(1, num_files);
    
    bead_parameter_fits = zeros(num_files, 6);
    bead_parameter_precisions = zeros(num_files, 6);
    bead_parameter_initial_guesses = zeros(num_files, 6);
    
    mask_lims_list = zeros(num_files, 2);
    
    applied_b_fields = zeros(num_files, 3);
    
    for i = 1:num_files
        
        file_name = list_of_files(i).name;
        
        load(file_name);
        
        measured_field = fit_data.measured_b_fields;
        
        fit_field_uT = fit_data.fit_bead_fields_uT;
        
        if i == 1
            
            measured_b_fields = zeros(num_files, ...
                size(measured_field, 1), size(measured_field, 2), 3);
            
            fit_b_fields = zeros(num_files, ...
                size(fit_field_uT, 1), size(fit_field_uT, 2), 3);
        end
        
        angles(i) = fit_data.angle;
        
        bead_parameter_fits(i, :) = fit_data.bead_parameter_fits;
        bead_parameter_precisions(i, :) = fit_data.bead_parameter_precisions;
        bead_parameter_initial_guesses(i, :) = fit_data.initial_guess;
        mask_lims_list(i, :) = fit_data.mask_lims;
        
        
        measured_b_fields(i, :, :, :) = fit_data.measured_b_fields;
        fit_b_fields(i, :, :, :) = fit_data.fit_bead_fields_uT;
        
        applied_b_fields(i, :) = fit_data.applied_field_uT;
    end
    
    [sorted_angles, sort_order] = sort(angles);
    
    angles = sorted_angles;
    
    temp = bead_parameter_fits(sort_order, :);
    bead_parameter_fits = temp;
    clearvars('temp');
    
    temp = bead_parameter_precisions(sort_order, :);
    bead_parameter_precisions = temp;
    clearvars('temp');
    
    temp = bead_parameter_initial_guesses(sort_order, :);
    bead_parameter_initial_guesses = temp;
    clearvars('temp');
    
    temp = mask_lims_list(sort_order, :);
    mask_lims_list = temp;
    clearvars('temp');
    
    temp = applied_b_fields(sort_order, :);
    applied_b_fields = temp;
    clearvars('temp');
    
    data.raw_data_file = fit_data.raw_data_file;
    data.angles = angles;
    
    data.fit_units = fit_data.fit_units;
    data.bead_parameter_fits = bead_parameter_fits;
    data.bead_parameter_precisions = bead_parameter_precisions;
    
    data.bead_parameter_initial_guesses = bead_parameter_initial_guesses;
    data.mask_lims_list = mask_lims_list;
    
    data.measure_b_fields_uT = measured_b_fields;
    data.fit_bead_fields_uT = fit_b_fields;
    
    data.applied_fields_uT = applied_b_fields;
    
    
    data_save_file_name = full_sweepfit_filename;
    
    cd(current_folder_full_path)    
    
    save(data_save_file_name, 'data');
    
else
    load(full_sweepfit_filename)
end

%% plot bead moment and applied field over rotation


% pulling out applied fields

num_fits = length(data.angles);

xy_projection_applied_field_dirs = zeros(num_fits, 2);

for i = 1:num_fits
    
    angle = data.angles(i);
    original_index = find(round(data.angle_space) == round(angle));
    
    applied_field = data.mean_fields_uT(original_index, :);
    
    xy_projection = applied_field(1:2);
    
    normed_xy_projection = xy_projection / norm(xy_projection);
    
    xy_projection_applied_field_dirs(i, :) = normed_xy_projection;
    
end

% pulling out moments

fit_phis = data.bead_parameter_fits(:, 5);

fit_moment_signs = sign(data.bead_parameter_fits(:, 6));

fit_moment_signs = ones(size(fit_moment_signs));

moment_x_components = cos(fit_phis).*fit_moment_signs;
moment_y_components = sin(fit_phis).*fit_moment_signs;

figure(6);
clf
for j = 1:num_fits
    
    moment_vector_length_rescale = 3/4;
    applied_field_vector_length_rescale = 2/3;
    
    hold on
    quiver( j, 0, ...
        moment_vector_length_rescale*moment_x_components(j), ...
        moment_vector_length_rescale*moment_y_components(j), '-b')
    quiver( j, 0, ...
        applied_field_vector_length_rescale*xy_projection_applied_field_dirs(j,1), ...
        applied_field_vector_length_rescale*xy_projection_applied_field_dirs(j,2), '-r')
    
end

xticks([])
yticks([])
axis([0 (num_fits+1) -1 1])
pbaspect([(num_fits+1) 2 1])
xticks(1:num_fits)
set(gca,'TickLength',[0, 0])
xticklabels(arrayfun(@num2str, data.angles', 'UniformOutput', 0))
legend({'magnetic moment projection', 'applied magnetic field projection'}, 'Location', 'southoutside')
    
%% plot rotation sequence (just fields)

num_fits = length(data.angles);

fit_sidelength = size(data.fit_bead_fields_uT, 2);


combined_fit_figure_x = zeros(fit_sidelength,fit_sidelength*num_fits);
combined_data_figure_x = zeros(fit_sidelength,fit_sidelength*num_fits);


combined_fit_figure_y = zeros(fit_sidelength,fit_sidelength*num_fits);
combined_data_figure_y = zeros(fit_sidelength,fit_sidelength*num_fits);

for j = 1:num_fits
    
    to_plot_fit_x = squeeze(data.fit_bead_fields_uT(j, :, :, 1));
    to_plot_measured_x = squeeze(data.measure_b_fields_uT(j, :, :, 1));
    
    
    to_plot_fit_y = squeeze(data.fit_bead_fields_uT(j, :, :, 2));
    to_plot_measured_y = squeeze(data.measure_b_fields_uT(j, :, :, 2));
    
    
    combined_fit_figure_x(:, 1+(j-1)*fit_sidelength:fit_sidelength*j) = to_plot_fit_x;
    combined_data_figure_x(:, 1+(j-1)*fit_sidelength:fit_sidelength*j) = to_plot_measured_x;
    
    combined_fit_figure_y(:, 1+(j-1)*fit_sidelength:fit_sidelength*j) = to_plot_fit_y;
    combined_data_figure_y(:, 1+(j-1)*fit_sidelength:fit_sidelength*j) = to_plot_measured_y;
    
      
end


plotting_x = false;

if plotting_x
    
    component_string = 'Bx'; 
    
    combined_fit_to_plot = combined_fit_figure_x;
    combined_measured_to_plot = combined_data_figure_x;
    
else
    
    component_string = 'By';
    
    combined_fit_to_plot = combined_fit_figure_y;
    combined_measured_to_plot = combined_data_figure_y;
    
end

caxis_range = [-2 2];


figure(4);
clf;

subplot(3,1,1)
imagesc(combined_fit_to_plot)
pbaspect([ fit_sidelength*num_fits fit_sidelength 1 ])
colormap(linspecer)
c = colorbar;
ylabel(c, [component_string ' (\muT)']);
caxis(caxis_range)
xticks(1+floor(fit_sidelength/2):fit_sidelength:floor((num_fits+1/2)*fit_sidelength))
set(gca,'TickLength',[0, 0])
xticklabels(arrayfun(@num2str, data.angles', 'UniformOutput', 0))
xlabel('MT azimuthal angle')
title(['Fit ' component_string ' rotation sequence'])
yticks([])
rectangle('Position', [5, 5, 10, 1], 'EdgeColor', 'k', 'FaceColor', 'k')

subplot(3,1,2)


imagesc(combined_measured_to_plot)

pbaspect([ fit_sidelength*num_fits fit_sidelength 1 ])
colormap(linspecer)
c = colorbar;
ylabel(c, [component_string ' (\muT)']);
caxis(caxis_range)

xticks(1+floor(fit_sidelength/2):fit_sidelength:floor((num_fits+1/2)*fit_sidelength))

set(gca,'TickLength',[0, 0])
xticklabels(arrayfun(@num2str, data.angles', 'UniformOutput', 0))
xlabel('MT azimuthal angle')
title([component_string ' rotation sequence'])
yticks([])
rectangle('Position', [5, 5, 10, 1], 'EdgeColor', 'k', 'FaceColor', 'k')

subplot(3,1,3)

imagesc(medfilt2(combined_measured_to_plot))

pbaspect([ fit_sidelength*num_fits fit_sidelength 1 ])
colormap(linspecer)
c = colorbar;
ylabel(c, 'By (\muT)');
caxis(caxis_range)

xticks(1+floor(fit_sidelength/2):fit_sidelength:floor((num_fits+1/2)*fit_sidelength))

set(gca,'TickLength',[0, 0])
xticklabels(arrayfun(@num2str, data.angles', 'UniformOutput', 0))
xlabel('MT azimuthal angle')
title([component_string ' rotation sequence (filtered)'])
yticks([])
rectangle('Position', [5, 5, 10, 1], 'EdgeColor', 'k', 'FaceColor', 'k')

%% plot rotation sequence (with overlaid projected bead moment and applied field)


% pulling out applied fields

num_fits = length(data.angles);

xy_projection_applied_field_dirs = zeros(num_fits, 2);

for i = 1:num_fits
    
    angle = data.angles(i);
    original_index = find(round(data.angle_space) == round(angle));
    
    applied_field = data.mean_fields_uT(original_index, :);
    
    xy_projection = applied_field(1:2);
    
    normed_xy_projection = xy_projection / norm(xy_projection);
    
    xy_projection_applied_field_dirs(i, :) = normed_xy_projection;
    
end

% pulling out moments

fit_phis = data.bead_parameter_fits(:, 5);

fit_moment_signs = sign(data.bead_parameter_fits(:, 6));

fit_moment_signs = ones(size(fit_moment_signs));

moment_x_components = cos(fit_phis).*fit_moment_signs;
moment_y_components = sin(fit_phis).*fit_moment_signs;

num_fits = length(data.angles);

fit_sidelength = size(data.fit_bead_fields_uT, 2);


combined_fit_figure_x = zeros(fit_sidelength,fit_sidelength*num_fits);
combined_data_figure_x = zeros(fit_sidelength,fit_sidelength*num_fits);


combined_fit_figure_y = zeros(fit_sidelength,fit_sidelength*num_fits);
combined_data_figure_y = zeros(fit_sidelength,fit_sidelength*num_fits);

for j = 1:num_fits
    
    to_plot_fit_x = squeeze(data.fit_bead_fields_uT(j, :, :, 1));
    to_plot_measured_x = squeeze(data.measure_b_fields_uT(j, :, :, 1));
    
    
    to_plot_fit_y = squeeze(data.fit_bead_fields_uT(j, :, :, 2));
    to_plot_measured_y = squeeze(data.measure_b_fields_uT(j, :, :, 2));
    
    
    combined_fit_figure_x(:, 1+(j-1)*fit_sidelength:fit_sidelength*j) = to_plot_fit_x;
    combined_data_figure_x(:, 1+(j-1)*fit_sidelength:fit_sidelength*j) = to_plot_measured_x;
    
    combined_fit_figure_y(:, 1+(j-1)*fit_sidelength:fit_sidelength*j) = to_plot_fit_y;
    combined_data_figure_y(:, 1+(j-1)*fit_sidelength:fit_sidelength*j) = to_plot_measured_y;
    
      
end


plotting_x = false;

if plotting_x
    
    component_string = 'Bx'; 
    
    combined_fit_to_plot = combined_fit_figure_x;
    combined_measured_to_plot = combined_data_figure_x;
    
else
    
    component_string = 'By';
    
    combined_fit_to_plot = combined_fit_figure_y;
    combined_measured_to_plot = combined_data_figure_y;
    
end

caxis_range = [-2 2];


figure(3);
clf;

subplot(3,1,1)
imagesc(combined_fit_to_plot)
pbaspect([ fit_sidelength*num_fits fit_sidelength 1 ])
colormap(linspecer)
c = colorbar;
ylabel(c, [component_string ' (\muT)']);
caxis(caxis_range)
xticks(1+floor(fit_sidelength/2):fit_sidelength:floor((num_fits+1/2)*fit_sidelength))
set(gca,'TickLength',[0, 0])
xticklabels(arrayfun(@num2str, data.angles', 'UniformOutput', 0))
xlabel('MT azimuthal angle')
title(['Fit ' component_string ' rotation sequence'])
yticks([])
rectangle('Position', [5, 5, 10, 1], 'EdgeColor', 'k', 'FaceColor', 'k')

subplot(3,1,2)

imagesc(combined_measured_to_plot)

pbaspect([ fit_sidelength*num_fits fit_sidelength 1 ])
colormap(linspecer)
c = colorbar;
ylabel(c, [component_string ' (\muT)']);
caxis(caxis_range)

xticks(1+floor(fit_sidelength/2):fit_sidelength:floor((num_fits+1/2)*fit_sidelength))

set(gca,'TickLength',[0, 0])
xticklabels(arrayfun(@num2str, data.angles', 'UniformOutput', 0))
xlabel('MT azimuthal angle')
title([component_string ' rotation sequence'])
yticks([])
rectangle('Position', [5, 5, 10, 1], 'EdgeColor', 'k', 'FaceColor', 'k')

subplot(3,1,3)

imagesc(medfilt2(combined_measured_to_plot))

pbaspect([ fit_sidelength*num_fits fit_sidelength 1 ])
colormap(linspecer)
c = colorbar;
ylabel(c, 'By (\muT)');
caxis(caxis_range)

xticks(1+floor(fit_sidelength/2):fit_sidelength:floor((num_fits+1/2)*fit_sidelength))

set(gca,'TickLength',[0, 0])
xticklabels(arrayfun(@num2str, data.angles', 'UniformOutput', 0))
xlabel('MT azimuthal angle')
title([component_string ' rotation sequence (filtered)'])
yticks([])
rectangle('Position', [5, 5, 10, 1], 'EdgeColor', 'k', 'FaceColor', 'k')

frame_centers = 1+floor(fit_sidelength/2):fit_sidelength:floor((num_fits+1/2)*fit_sidelength);

for k = 1:3
    subplot(3,1,k)
    for j = 1:num_fits
        
        % negative sign is because quiver and imagesc have inverted y
        % coordinates
        moment_vector_length_rescale = -5;
        applied_field_vector_length_rescale = -3;
        
        x_val = frame_centers(j)+floor(fit_sidelength/4);
        y_val = fit_sidelength - abs(moment_vector_length_rescale);
        
        hold on
        quiver( x_val, y_val, ...
            moment_vector_length_rescale*moment_x_components(j), ...
            moment_vector_length_rescale*moment_y_components(j), '-b', ...
            'LineWidth', 1)
        quiver( x_val, y_val, ...
            applied_field_vector_length_rescale*xy_projection_applied_field_dirs(j,1), ...
            applied_field_vector_length_rescale*xy_projection_applied_field_dirs(j,2), '-r', ...
            'LineWidth', 1)
        
    end
    
    legend({'magnetic moment projection', 'applied magnetic field projection'}, 'Location', 'southoutside')
end

%% saving the figures that were generated

