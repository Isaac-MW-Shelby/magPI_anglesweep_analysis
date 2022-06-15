function [fit_outputs, fit_precisions] = fit_applied_field(angle_data, ... 
                        shift_magnitudes, orientations, rescale_vector, debug)
                    
[~, max_ind] = max(shift_magnitudes(:, 4));

shift_vals_at_max = shift_magnitudes(max_ind, :);

theta_val_of_max = angle_data(max_ind);

B0_guess_inner_two = (shift_vals_at_max(1)+shift_vals_at_max(2))*(3/4);
B0_guess_outer_two = (shift_vals_at_max(4) - shift_vals_at_max(3))*3/4;


B_split_mag_guess = mean([B0_guess_inner_two B0_guess_outer_two]); 
B_split_component_guess = B_split_mag_guess/sqrt(3);

B_tweezer_mag_guess = (shift_vals_at_max(4) - 3*shift_vals_at_max(3))*((-1/4)*(sqrt(3/2)));

B_tweezer_phi_offset_guess = mod(pi/4*180/pi - theta_val_of_max, 360);

% T0_manual = [1.4*10^4, 1.7*10^4, 0.8*10^4, 3.5*10^4, 131, 3.7*10^3, 131];
% 
% T0 = T0_manual;  

T0 = [B_split_component_guess B_split_component_guess ...
      B_split_component_guess ...
      B_tweezer_mag_guess B_tweezer_phi_offset_guess ...
      B_tweezer_mag_guess/10 B_tweezer_phi_offset_guess];
  


lb = zeros(size(T0));
lb(1) = 0;
lb(2) = 0;
lb(3) = 0;
lb(4) = 0;
lb(5) = -720;
lb(6) = 0;
lb(7) = -720;



ub = zeros(size(T0));

ub(1) = 1*28000; % KHz
ub(2) = 1*28000; % KHz
ub(3) = 1*28000; % KHz
ub(4) = 3*28000; % KHz
ub(5) = 720;
ub(6) = 1*28000; % KHz
ub(7) = 720;

T0 = T0 ./ rescale_vector;
lb = lb ./ rescale_vector;
ub = ub ./ rescale_vector;

options.OptimalityTolerance = 1*10^-10;
options.FunctionTolerance = 1*10^-10;
options.Display = 'off';
    
[Tfit,~,~,~,~,~,jacob] = lsqnonlin( @tmpFnt , T0, lb, ub,options );

fit_outputs = Tfit.*rescale_vector;

fisher_matrix = jacob'*jacob;
correlation_matrix = inv(fisher_matrix);
unscaled_precisions = (diag(correlation_matrix)).^(1/2);

fit_precisions = unscaled_precisions(1:length(rescale_vector)).*rescale_vector;

    function dy = tmpFnt(T)
        
        dy = zeros(size(shift_magnitudes));
        
        for i = 1:length(angle_data)
            
            magnetic_field = applied_magnetic_field(T.*rescale_vector, ...
                 angle_data(i));
            
            shifts = spectrum_from_b_field(magnetic_field,orientations);
            
            dy(i, :) = sort(abs(shift_magnitudes(i, :))) - sort(abs(shifts))';
        end
        
        dy = dy(:);
        
        if debug
            
            applied_field_parameters = T.*rescale_vector;
            
            num_guess_angles = 100;
            
            theory_angle_space = linspace(0, 360, num_guess_angles);
            
            field = zeros(length(theory_angle_space), 3);
            
            for i = 1:num_guess_angles
                field(i, :) = ...
                    applied_magnetic_field(applied_field_parameters, ...
                    theory_angle_space(i));
            end
            
            
            figure(1011)
            clf
            for i = 1:4
                plot( angle_data, shift_magnitudes(:, i), '+','LineWidth',1, 'Color','b');
                hold on
            end
            title('resonance shifts vs tweezer magnetic field angle with overlaid initial guess')
            xlabel('tweezer magnetic field angle')
            ylabel('resonance shift (kHz)')
            xlim([-5 365])
            xticks(angle_data)
            
            
            guess_applied_B_field = zeros(num_guess_angles, 3);
            guess_applied_B_field(:, 1) =  field(:, 1);
            guess_applied_B_field(:, 2) =  field(:, 2);
            guess_applied_B_field(:, 3) =  field(:, 3);
            

            
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
            
            figure(1010001)
            clf
            subplot(1,3,1)
            plot(guess_applied_B_field(:, 1))
            subplot(1,3,2)
            plot(guess_applied_B_field(:, 2))
            subplot(1,3,3)
            plot(guess_applied_B_field(:, 3))
            
            pause
            
        end
        
    end


                    
end

                    
            