

% run this in the folder with the B_field_maps.mat file 


%% fitting a single angle


angle_sweep_location = 'B_field_maps.mat';

angle_index_input = 5;


[bead_parameters] = ...
    fitting_bead_parameters_from_b_fields_wrapper(angle_sweep_location, angle_index_input);


%% fitting the whole sweep 


angle_sweep_location = 'B_field_maps.mat';

load(angle_sweep_location, 'data')

num_angles = length(data.angle_space);

angles_indices_fitting = 1:num_angles;

% can manually change angle_indices_fitting to fit specific angles
% angles_indices_fitting = [];

for angle_index = angles_indices_fitting
    
    disp(['starting fitting bead parameters for \phi_{MT} = ' ...
        num2str(data.angle_space(angle_index))])
    [bead_parameters] = ...
    fitting_bead_parameters_from_b_fields_wrapper(angle_sweep_location, angle_index);
    disp(['done fitting bead parameters for \phi_{MT} = ' ...
        num2str(data.angle_space(angle_index))])
    
end

