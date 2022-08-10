%% initial flags / inputs

clear all %#ok<CLALL>
close all

simulated_data = false;

do_save = true;

debug = true;

work_computer = true;

convert_sensitivities_to_precisions = true;

turbobead_data = true;

if work_computer
    slash = '/';
else 
    slash = '\';
end

loading_files_from_current_folder = true;

debug_figure_tracker = 1000;

%NV orientations
normalization = 1/(sqrt(3));

o1 = normalization * [1 1 1];
o2 = normalization * [-1 -1 1];
o3 = normalization * [-1 1 -1];
o4 = normalization * [1 -1 -1];
orientations = [o1; o2 ;o3 ; o4];

% gyromagnetic ratios
kHz_per_mT = 28*10^3;
kHz_per_uT = 28;

% for gating bad fits
max_shift_size = 0.25; % GHz corresponds to ~8.9 mT

figure_number_tracker = 1;

%% select data files to load

if simulated_data 
    
    prefix = '/home/ishelby/Documents/MATLAB/magnetic_torque_project/all_four_NV_vector_sensing/simulated_data';
    experimental_data_file_name = 'simulated_torquing_data.mat';
    %    experimental_data_file_name = 'simulated_torquing_data_nonzero_tilt.mat';
    
    load([prefix '/' experimental_data_file_name])
    
    shift_data = data.resonance_splittings;
    
    size_shift_data = size(shift_data);
    
    angle_space = linspace(0, 360, size_shift_data(1) + 1);
    
    angle_space = angle_space(1:end-1);
    
else
    
    angle_sweep_prefix = 'anglesweep_WFfullresonancemaps_';
    
    if ~loading_files_from_current_folder
    
    
        %experimental_data_file_name = ...
        %                           'anglesweep_WFfullresonancemaps.mat';
        experimental_data_file_name = '020822_E1421_r1';

        
        if work_computer
            prefix = '/home/ishelby/Documents/MATLAB/magnetic_torque_project/all_four_NV_vector_sensing/real_data'; %#ok<*UNRCH>
            data_location = [prefix '/files_for_analysis_' experimental_data_file_name];
            cd(data_location);
        else
            
            prefix = 'C:\Users\Isaac\Documents\MATLAB\dna-torquing-project\all_four_NV_vector_sensing\real_data';
            load([prefix '\' experimental_data_file_name])
        end
        
        load([angle_sweep_prefix experimental_data_file_name '.mat'])
    
    else
        
        load_prefix = [angle_sweep_prefix '*'];
        
        main_file = dir(load_prefix);
        
        load(main_file.name)
        
        experimental_data_file_name = main_file.name(length(angle_sweep_prefix)+1:end-4);
        
        
    end
    
end

%% pulling data from loaded file
    
if exist('croppedmaps', 'var')
    resonance_locations = croppedmaps;
    disp('assigned croppedmaps to resonance_locations')
elseif exist('anglesweep_WFfullresonancemaps', 'var')
    resonance_locations =  anglesweep_WFfullresonancemaps;
    disp('assigned anglesweep_WFfullresonancemaps to resonance_locations')
end


if resonance_locations(1,1,1,1) > resonance_locations(1,1,1,8)
    resonance_locations = flip(resonance_locations,4);
end

size_data = size(resonance_locations);

shift_data = zeros(size_data(1), size_data(2), size_data(3), size_data(4)/2);

for i = 1:4
    
    shift_data(:, :, :, i) = resonance_locations(:, :, :, 9 - i) - resonance_locations(:, :, :, i);
    
end


angle_indices_used = 1:size(shift_data,1);


temp = shift_data;
shift_data = shift_data(angle_indices_used, :, :, :);

size_shift_data = size(shift_data);

for i = 1:size_shift_data(1)
    for j = 1:size_shift_data(2)
        for k = 1: size_shift_data(3)
            shift_data(i, j, k, :) = sort(shift_data(i, j, k, :));
        end
    end
end



if exist('angles', 'var')
    angle_space = angles(angle_indices_used);
else
    angle_space = linspace(0, 360, size_shift_data(1));
end

%% loading in matrices outputs for sensitivities to give precisions


if convert_sensitivities_to_precisions

    extra_matrix_directory = dir('*magangle*');
    
    num_angles = size_shift_data(1);
    
    precisions_uT = zeros([size_shift_data(1:3) 3]);
    
    pulled_angles = zeros(1, size_shift_data(1));
    
    
    for i = 1:num_angles
        
        matrix_filename = extra_matrix_directory(i).name;
        load(matrix_filename);
        
        angle_of_mat_str = extractAfter(matrix_filename, 'magangle');
        angle_of_mat = round(str2num(angle_of_mat_str(1:end-4))); %#ok<ST2NM>
        
        
        
        matrix_angle_index = find(angle_of_mat == round(angle_space));
        
        pulled_angles(matrix_angle_index) = angle_of_mat;
        
        conversion_factor = (10^6)*0.55/0.6679;
        % converts GHz/sqrt(Hz) (sens) -> kHz/sqrt(Hz) (sens) -> uT (precision)
        % .55 uT for .6679 kHz/sqrt(Hz) found phemonenologically from empty FOV
        
        for j = 1:3
            
            precisions_uT(i, :, :, j) = conversion_factor*c.sensitivities;
            
        end
    end

else
    
    precisions_uT = 0.6*ones([size_shift_data(1:3) 3]); % 0.6 uT is from experiment
    
end

%% generates parent folder and moves data to that

current_location = pwd;

location_of_slashes = find(current_location == '/');

last_slash_location = location_of_slashes(end);
penulimate_slash_location = location_of_slashes(end-1);

current_folder_name = current_location(last_slash_location+1:end);

parent_folder_full_path = current_location(1:last_slash_location-1);

parent_folder_name = current_location(penulimate_slash_location+1:last_slash_location-1);

cd(parent_folder_full_path);

data_file_folder_prefix = 'files_for_analysis_';

new_parent_folder_name = current_folder_name(length(data_file_folder_prefix)+1:end);

if length(parent_folder_name) == length(new_parent_folder_name)

    already_run = parent_folder_name == new_parent_folder_name;
else
    already_run = false;
    
end


if ~already_run

    if ~exist(new_parent_folder_name, 'dir')
        mkdir(new_parent_folder_name)
    end
    
    movefile(current_folder_name, new_parent_folder_name)
    
    cd(new_parent_folder_name)
end

%% catching weirdness in the data (usually superfluous now)


negative_shifts = shift_data < 0;

negative_shift_locations = find(negative_shifts == 1);

[neg_shift_angle_inds, neg_shift_x_pixel_inds, neg_shift_y_pixel_inds, ...
    neg_shift_order_inds] = ...
    ind2sub(size_shift_data, negative_shift_locations);

if simulated_data
    max_shift_size = max_shift_size*10^6;
    
end

overly_large_shifts = abs(shift_data) > max_shift_size;

overly_large_shift_locations = find(overly_large_shifts == 1);

[ol_shift_angle_inds, ol_shift_x_pixel_inds, ol_shift_y_pixel_inds, ...
    ol_shift_order_inds] = ...
    ind2sub(size_shift_data, overly_large_shift_locations);

data_to_use = (~negative_shifts) .* (~overly_large_shifts);

data_to_throw_out = ~data_to_use;

data_to_throw_out_locations = find(data_to_throw_out == 1);

[NaN_angle_inds, NaN_x_pixel_inds, NaN_y_pixel_inds, ...
    NaN_order_inds] = ...
    ind2sub(size_shift_data, data_to_throw_out_locations);

shift_data_NaNed_bad_data = shift_data;

for i = 1:length(NaN_angle_inds)
    
    angle = NaN_angle_inds(i);
    x_pixel = NaN_x_pixel_inds(i);
    y_pixel = NaN_y_pixel_inds(i);
    
    shift_data_NaNed_bad_data(angle, x_pixel, y_pixel, :) = [NaN NaN NaN NaN];
end

unused_data = isnan(shift_data_NaNed_bad_data);

fraction_of_usable_data = 1-sum(unused_data(:)) / numel(shift_data_NaNed_bad_data);

if fraction_of_usable_data < 1
    
    disp(['threw out negative/overly large shifts. kept ' num2str(fraction_of_usable_data*100) ' percent'])
else
    disp('no negative/overly large shifts')
end

shift_data_to_analyze = shift_data_NaNed_bad_data;

%% calculate mean shifts and plot vs tweezer angle


% first we find the average shift for all 4 orientations at each angle
% note that at this point we do not which orientation is, and we don't know their signs

if ~simulated_data
    shift_data_to_analyze = shift_data_to_analyze*10^6; % converts GHz to kHz
end

mean_shifts = zeros(size_shift_data(1), 4);

for i = 1:size_shift_data(1)
    
    mean_resonance_shifts = zeros(4, 1);
    
    for j = 1:4
        
        data_to_mean = squeeze(shift_data_to_analyze(i, :, :, j));
        mean_resonance_shifts(j) =  nanmean(data_to_mean(:));
        
    end
    
    
    mean_shifts(i, :) = sort(mean_resonance_shifts);
    
    
end


current_fig = figure(figure_number_tracker);
figure_number_tracker = figure_number_tracker + 1;

clf

for i = 1:4
    plot( angle_space, mean_shifts(:, i), '+','LineWidth',1, 'Color','b');
    hold on
end
hold off
fig_title = 'mean resonance shifts vs tweezer magnetic field angle';
current_fig.Name = fig_title;
title(fig_title)
xlabel('tweezer magnetic field angle (deg)')
ylabel('resonance shift (kHz)')
xlim([-5 365])
xticks(angle_space)

%% get initial guess for applied B fields


[max_val, max_ind] = max(mean_shifts(:, 4));

shift_vals_at_max = mean_shifts(max_ind, :);

theta_val_of_max = angle_space(max_ind);

B0_guess_inner_two = (shift_vals_at_max(1)+shift_vals_at_max(2))*(3/4);
B0_guess_outer_two = (shift_vals_at_max(4) - shift_vals_at_max(3))*3/4;
B0_guess = mean([B0_guess_inner_two B0_guess_outer_two]);

B1_guess = (shift_vals_at_max(4) - 3*shift_vals_at_max(3))*((-1/4)*(sqrt(3/2)));

offset_guess = pi/4*180/pi - theta_val_of_max;

tilt_guess = pi/2;

%% fit full angle sweep for resonance assignment

[resonance_assignment, fit_outputs, fit_precisions] = ...
    fit_resonance_assignment(angle_space, mean_shifts, ...
    orientations);

% check strings are hard coded fixes to the resonance assignment

check_string = 'anglesweep_WFfullresonancemaps_UnconstrainedST_010422_r31.mat';
length_restricted_check_string = check_string(1:length(experimental_data_file_name));

if experimental_data_file_name == length_restricted_check_string
    
    resonance_assignment(3,1) = abs(resonance_assignment(3,1));
    resonance_assignment(6,1) = abs(resonance_assignment(6,1));
    
end


check_string = '020822_E1421_r1';
length_restricted_check_string = check_string(1:length(experimental_data_file_name));

if experimental_data_file_name == length_restricted_check_string
    
    temp = resonance_assignment(1,2);
    resonance_assignment(1,2) = resonance_assignment(1,3);
    resonance_assignment(1,3) = temp;
    
    temp = resonance_assignment(end,2);
    resonance_assignment(end,2) = resonance_assignment(end,3);
    resonance_assignment(end,3) = temp;
    
end

B_split_x = fit_outputs(1);
B_split_y = fit_outputs(2);
B_split_z = fit_outputs(3);
splitting_field = [B_split_x B_split_y B_split_z];

B_mt_magnitude_in_plane = fit_outputs(4);
B_mt_phase_in_plane = fit_outputs(5);
B_mt_magnitude_out_of_plane = fit_outputs(6);
B_mt_phase_out_of_plane = fit_outputs(7);

splitting_field_uT = [B_split_x B_split_y B_split_z] / kHz_per_uT;
B_mt_magnitude_in_plane_uT = B_mt_magnitude_in_plane / kHz_per_uT;
B_mt_magnitude_out_of_plane_uT = B_mt_magnitude_out_of_plane / kHz_per_uT;

%% plot splittings due to fit applied field against real splittings

num_guess_angles = 201;

theory_angle_space = linspace(0, 360, num_guess_angles);

Bx = splitting_field(1) + B_mt_magnitude_in_plane*cos((pi/180)*(theory_angle_space + B_mt_phase_in_plane));
By = splitting_field(2) + B_mt_magnitude_in_plane*sin((pi/180)*(theory_angle_space + B_mt_phase_in_plane));
Bz = splitting_field(3) + B_mt_magnitude_out_of_plane*sin((pi/180)*(theory_angle_space + B_mt_phase_out_of_plane));

current_fig = figure(figure_number_tracker);
figure_number_tracker = figure_number_tracker + 1;
for i = 1:4
    plot( angle_space, mean_shifts(:, i), '+','LineWidth',1, 'Color','b');
    hold on
end
fig_title = 'fit mean applied field shifts and raw data';
current_fig.Name = fig_title;
title(fig_title)
xlabel('tweezer magnetic field angle (deg)')
ylabel('resonance shift (kHz)')
xlim([-5 365])
xticks(angle_space)


guess_applied_B_field = zeros(num_guess_angles, 3);
guess_applied_B_field(:, 1) =  Bx;
guess_applied_B_field(:, 2) =  By;
guess_applied_B_field(:, 3) =  Bz;

shifts = zeros(num_guess_angles, 4);


for t = 1:num_guess_angles
    
    for j = 1:4
        
        shifts(t, j) = 2*abs(dot(orientations(j, :), guess_applied_B_field(t, :)));
        
    end
    
end



o1_guess = plot(theory_angle_space, shifts(:, 1));
o2_guess = plot(theory_angle_space, shifts(:, 2));
o3_guess = plot(theory_angle_space, shifts(:, 3));
o4_guess = plot(theory_angle_space, shifts(:, 4));

legend([o1_guess, o2_guess, o3_guess, o4_guess], 'o1 guess', 'o2 guess', 'o3 guess', 'o4 guess')

hold off

%% generate B field maps

[per_pixel_B_fields, per_pixel_B_precisions] = ...
    generateBFieldMaps(shift_data*10^6, resonance_assignment);

mean_fields = zeros(size(per_pixel_B_fields,1), 3);

for i = 1:size(per_pixel_B_fields,1)
    
    Bx = per_pixel_B_fields(i, :, :, 1);
    mean_fields(i, 1) = mean(Bx(:));
    
    By = per_pixel_B_fields(i, :, :, 2);
    mean_fields(i, 2) = mean(By(:));
    
    Bz = per_pixel_B_fields(i, :, :, 3);
    mean_fields(i, 3) = mean(Bz(:));
end

mean_fields_uT = mean_fields / kHz_per_uT;

current_fig = figure(figure_number_tracker);
figure_number_tracker = figure_number_tracker + 1;

current_fig.Name = 'mean field components vs tweezer angle';
subplot(2,2,1)
plot(angle_space, mean_fields_uT(:, 1))
title('mean Bx')
xlabel('MT angle (degrees)')
ylabel('Bx (uT)')


subplot(2,2,2)
plot(angle_space, mean_fields_uT(:, 2))
title('mean By')
xlabel('MT angle (degrees)')
ylabel('By (uT)')


subplot(2,2,3)
plot(angle_space, mean_fields_uT(:, 3))
title('mean Bz')
xlabel('MT angle (degrees)')
ylabel('Bz (uT)')


subplot(2,2,4)
plot(angle_space, (mean_fields_uT(:, 1).^2 + ...
                    mean_fields_uT(:, 2).^2 + ...
                        mean_fields_uT(:, 3).^2).^(1/2));
title('magnitude of mean field')
xlabel('MT angle (degrees)')
ylabel('B (uT)')

%% generate mean subtracted bead field maps

mean_subtracted_fields_uT = zeros(size(per_pixel_B_fields));

for i = 1:size(per_pixel_B_fields,1)
    
    
    Bx = (squeeze(per_pixel_B_fields(i, :, :, 1))-mean_fields(i,1))/ kHz_per_uT;
    By = (squeeze(per_pixel_B_fields(i, :, :, 2))-mean_fields(i,2))/ kHz_per_uT;
    Bz = (squeeze(per_pixel_B_fields(i, :, :, 3))-mean_fields(i,3))/ kHz_per_uT;
    
    mean_subtracted_fields_uT(i, :, :, 1) = Bx;
    mean_subtracted_fields_uT(i, :, :, 2) = By;
    mean_subtracted_fields_uT(i, :, :, 3) = Bz;
    
    magnitude = (Bx.^2 + By.^2 + Bz.^2).^(1/2);
    
    
    if turbobead_data
        caxis_range = [-56 56]/ kHz_per_uT; % converted to uT
        magnitude_caxis_range = [0 100] / kHz_per_uT; % converted to uT
    else
        caxis_range = [-200 200]/ kHz_per_uT; % converted to uT
        magnitude_caxis_range = [0 300] / kHz_per_uT; % converted to uT
    end
    
    
    
    current_fig = figure(figure_number_tracker);
    figure_number_tracker = figure_number_tracker + 1;
    
    
    clf
    
    current_fig.Name = ['bead field components MT angle ' num2str(angle_space(i))];
    
    subplot(2,2,1)
    imagesc(Bx)
    pbaspect([1 1 1])
    title(['Bx \phi_{MT} = ' num2str(angle_space(i))])
    c = colorbar;
    c.Label.String = 'Bx (uT)';
    colormap(linspecer)
    caxis(caxis_range)
    xticks([])
    yticks([])
    
    subplot(2,2,2)
    imagesc(By)
    pbaspect([1 1 1])
    title(['By \phi_{MT} = ' num2str(angle_space(i))])
    c = colorbar;
    c.Label.String = 'By (uT)';
    colormap(linspecer)
    caxis(caxis_range)
    xticks([])
    yticks([])
    
    subplot(2,2,3)
    imagesc(Bz)
    pbaspect([1 1 1])
    title(['Bz \phi_{MT} = ' num2str(angle_space(i))])
    c = colorbar;
    c.Label.String = 'Bz (uT)';
    colormap(linspecer)
    caxis(caxis_range)
    xticks([])
    yticks([])
    
    subplot(2,2,4)
    imagesc(magnitude)
    pbaspect([1 1 1])
    title(['B magnitude \phi_{MT} = ' num2str(angle_space(i))])
    c = colorbar;
    c.Label.String = 'B (uT)';
    colormap(linspecer)
    caxis(magnitude_caxis_range)
    xticks([])
    yticks([])
    
    pause(1)
    
end

%% background distributions (for empty FOV's, commented out)


% mean_subtracted_field_mag_uT = zeros(size(mean_subtracted_fields_uT, 1:3));
%
% for i = 1:size(per_pixel_B_fields,1)
%
%    angle = angle_space(i);
%
%    Bx = mean_subtracted_fields_uT(i, :, :, 1);
%    By = mean_subtracted_fields_uT(i, :, :, 2);
%    Bz = mean_subtracted_fields_uT(i, :, :, 3);
%
%    mean_subtracted_field_mag_uT(i, :, :) = (Bx.^2 + By.^2 + Bz.^2).^(1/2);
%    figure(10000 +i)
%    clf
%    subplot(2,2,1)
%    histogram(Bx(:))
%    title(['Mean subtracted Bx distribution (uT) at angle ' num2str(angle)])
%    xlabel('Mean subtracted Bx (uT)')
%    ylabel('number of pixels')
%
%    subplot(2,2,2)
%    histogram(By(:))
%    title('Mean subtracted By distribution (uT)')
%    xlabel('Mean subtracted By (uT)')
%    ylabel('number of pixels')
%
%    subplot(2,2,3)
%    histogram(Bz(:))
%    title('Mean subtracted Bz distribution (uT)')
%    xlabel('Mean subtracted Bz (uT)')
%    ylabel('number of pixels')
%
%    subplot(2,2,4)
%    histogram(mean_subtracted_field_mag_uT(i, :, :))
%    title('Mean subtracted B magnitude distribution (uT)')
%    xlabel('Mean subtracted B magnitude (uT)')
%    ylabel('number of pixels')
%
%
% end
%
% start_index = 8;
%
%
% figure(11901)
% clf
% subplot(2,2,1)
% histogram(mean_subtracted_fields_uT(start_index:end, :, :, 1))
% title(['Mean subtracted Bx distribution (uT) at angle ' num2str(angle)])
% xlabel('Mean subtracted Bx (uT)')
% ylabel('number of pixels')
% xlim([-2 2])
% ylim([0 2500])
%
% subplot(2,2,2)
% histogram(mean_subtracted_fields_uT(start_index:end, :, :, 2))
% title('Mean subtracted By distribution (uT)')
% xlabel('Mean subtracted By (uT)')
% ylabel('number of pixels')
% xlim([-2 2])
% ylim([0 2500])
%
% subplot(2,2,3)
% histogram(mean_subtracted_fields_uT(start_index:end, :, :, 3))
% title('Mean subtracted Bz distribution (uT)')
% xlabel('Mean subtracted Bz (uT)')
% ylabel('number of pixels')
% xlim([-2 2])
% ylim([0 2500])
%
% subplot(2,2,4)
% histogram(mean_subtracted_field_mag_uT(start_index:end, :, :))
% title('Mean subtracted B magnitude distribution (uT)')
% xlabel('Mean subtracted B magnitude (uT)')
% ylabel('number of pixels')
% xlim([0 3])
% ylim([0 2500])

%% saving stuff

% debugging stuff, leave commented out
% experimental_data_file_name = experimental_data_file_name(1:15);
% experimental_data_file_name = [experimental_data_file_name '_test'];

current_location = pwd;

if do_save
    
%     
%     if work_computer
%         prefix = '/home/ishelby/Documents/MATLAB/magnetic_torque_project/all_four_NV_vector_sensing/fit_outputs';
%         
%         slash = '/';
%         
%         full_save_folder_path = [prefix slash experimental_data_file_name];
%     else
%         prefix = 'C:\Users\Isaac\Documents\MATLAB\magnetic_torque_project_9_9_2021\magnetic_torque_project\all_four_NV_vector_sensing\fit_outputs';
%         slash = '\';
%         full_save_folder_path = [prefix slash experimental_data_file_name];
%     end
%     
%     full_save_folder_path = erase(full_save_folder_path, '.mat');
%     
%     if ~exist(full_save_folder_path, 'dir')
%         mkdir(full_save_folder_path)
%     end
%     
%     cd(full_save_folder_path);
    
    data.raw_data_name = experimental_data_file_name;
    data.angle_space = angle_space;
    data.shift_data = shift_data;
    data.mean_subtracted_fields_uT = mean_subtracted_fields_uT;
    data.precisions_uT = precisions_uT;
    data.mean_fields_uT = mean_fields_uT;
    
    data.splitting_field_uT = splitting_field_uT;
    data.B_mt_magnitude_in_plane_uT = B_mt_magnitude_in_plane_uT;
    data.B_mt_phase_in_plane = B_mt_phase_in_plane;
    data.B_mt_magnitude_out_of_plane_uT = B_mt_magnitude_out_of_plane_uT;
    data.B_mt_phase_out_of_plane = B_mt_phase_out_of_plane;
    
    data.resonance_assignment = resonance_assignment;
    
    save('B_field_maps.mat', 'data')
    
    % saving the figures
    
    figure_folder_name = 'figures_from_convert_resonance_maps_to_B_fields';
    
    
    if ~exist(figure_folder_name, 'dir')
        mkdir(figure_folder_name)
    end
    
    cd(figure_folder_name)
    
    filetype_list = {'png', 'fig', 'eps'};
    
    
    for i = 1:length(filetype_list)
        
        filetype = filetype_list{i};
        
        if ~exist(filetype, 'dir')
            mkdir(filetype);
        end
        
    end
    
    cd(current_location);
    
    for i = 1:figure_number_tracker-1
        
        current_fig = figure(i);



        fig_title = current_fig.Name;
        fig_title = replace(fig_title, ' ', '_');
        
        fig_title = [experimental_data_file_name '_' fig_title]; %#ok<AGROW>
        
        for j = 1:3
            
            filetype = filetype_list{j};
        
            full_save_string = [figure_folder_name slash filetype slash fig_title];
            doPageFormat(2*[5,3]);
        
            full_save_string = replace(full_save_string, '.', 'point');
        
            saveas(current_fig, full_save_string, filetype);
        end
        
    end
    
end


