function [resonance_assignment, fit_outputs, fit_precisions] = ...
    fit_resonance_assignment(angle_data, shift_magnitudes, orientations)


angle_rescale = 360;

field_rescale = 2*28000;

rescale_vector = [field_rescale, field_rescale, field_rescale, ...
                    field_rescale, angle_rescale, ...
                        field_rescale, angle_rescale];
                    
debug  = false;

[fit_outputs, fit_precisions] = fit_applied_field(angle_data, ... 
                        shift_magnitudes, orientations, rescale_vector, ...
                        debug);

resonance_assignment = zeros(size(shift_magnitudes));

splitting_signs = zeros(size(resonance_assignment));

for i=1:size(resonance_assignment,1)
    
    
    magnetic_field = applied_magnetic_field(fit_outputs, ...
        angle_data(i));
    
    shifts = spectrum_from_b_field(magnetic_field,orientations);
    
   [~, resonance_assignment(i, :)] = sort(abs(shifts));
   
   unsorted_signs = zeros(1,4);
   
   for j = 1:4
       unsorted_signs(j) = sign(dot(orientations(j, :), magnetic_field));
   end
   
   splitting_signs(i, :) = unsorted_signs(resonance_assignment(i, :));
   
   
end

resonance_assignment = resonance_assignment .* splitting_signs;

end