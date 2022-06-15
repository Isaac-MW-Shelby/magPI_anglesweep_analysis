function [applied_magnetic_field] = ...
    applied_magnetic_field(applied_field_parameters, mt_angle)

B_split_x = applied_field_parameters(1);
B_split_y = applied_field_parameters(2);
B_split_z = applied_field_parameters(3);

B_mt_magnitude_in_plane     = applied_field_parameters(4);
B_mt_phase_in_plane         = applied_field_parameters(5);
B_mt_magnitude_out_of_plane = applied_field_parameters(6);
B_mt_phase_out_of_plane     = applied_field_parameters(7);

B_splitting = [B_split_x B_split_y B_split_z];

B_mt_x = ...
    B_mt_magnitude_in_plane.*cos((pi/180).*(mt_angle + B_mt_phase_in_plane));
B_mt_y = ...
    B_mt_magnitude_in_plane.*sin((pi/180).*(mt_angle + B_mt_phase_in_plane));
B_mt_z = ...
    B_mt_magnitude_out_of_plane.*...
                    sin((pi/180).*(mt_angle + B_mt_phase_out_of_plane));

B_mt = [B_mt_x B_mt_y B_mt_z];

applied_magnetic_field = B_splitting + B_mt;


end

