

% run this in the folder with the B_field_maps.mat file 


%% fitting a single angle


angle_sweep_location = 'B_field_maps.mat';

angle_index_input = 1;

bead_height = 600;

[bead_parameters] = ...
        fitting_bead_parameters_from_b_fields_wrapper(angle_sweep_location,...
        angle_index_input, bead_height);
    
%% sweeping over fixed z heights

angle_sweep_location = 'B_field_maps.mat';

bead_heights = 800;

for bead_height = bead_heights
    [bead_parameters] = ...
        fitting_bead_parameters_from_b_fields_wrapper(angle_sweep_location,...
        angle_index_input, bead_height);
end

%% sweeping over gradient masking vals

angle_sweep_location = 'B_field_maps.mat';

gradient_vals = 24; % in nT/nm

gradient_vals = gradient_vals / 1000; % converts to uT/nm

for grad_val = gradient_vals
    [bead_parameters] = ...
        fitting_bead_parameters_from_b_fields_wrapper(angle_sweep_location,...
        angle_index_input, grad_val);
end


%% fitting the whole sweep 


angle_sweep_location = 'B_field_maps.mat';

load(angle_sweep_location, 'data')

num_angles = length(data.angle_space);

angles_indices_fitting = 1:num_angles;

% can manually change angle_indices_fitting to fit specific angles
angles_indices_fitting = [1 3 4 6 7 8 9 10 11 12 13 14];

bead_height = 600; % nm height of turbobead stuck on surface

for angle_index = angles_indices_fitting
    
    disp(['starting fitting bead parameters for \phi_{MT} = ' ...
        num2str(data.angle_space(angle_index))])
    [bead_parameters] = ...
    fitting_bead_parameters_from_b_fields_wrapper(angle_sweep_location, angle_index, bead_height);
    disp(['done fitting bead parameters for \phi_{MT} = ' ...
        num2str(data.angle_space(angle_index))])
    
end

