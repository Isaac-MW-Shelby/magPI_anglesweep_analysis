

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
end

%% plot bead moment and applied field over rotation

%% plot rotation sequence (just fields)

%% plot rotation sequence (with overlaid projected bead moment and applied field)



