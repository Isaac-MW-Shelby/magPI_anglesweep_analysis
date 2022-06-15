function [bead_parameter_output, precisions] = fitting_bead_parameters_from_b_fields( ...
    X,Y, nvDepth, initial_guess, measured_fields, ...
    measured_field_precisions, rescale_vector, ...
    gradient_value, conversion_to_uT, rescale_data_to_fit, ...
    pixel_side_length, subdivisions_per_axis, ...
    front_binning_matrix, back_binning_matrix, fixed_z_val, mask_max_dist, mask_min_dist, debug)


    options.Display = 'off';
    options.FunctionTolerance = 1*10^-12;
    options.OptimalityTolerance = 1*10^-10;
    options.StepTolerance = 1*10^-12;
    options.MaxFunEvals = 6*100*10;
    options.MaxIter = 400*10;
    options.Algorithm = 'levenberg-marquardt';
    
    T0 = initial_guess ./ rescale_vector;
    
    % assumes cropped well enough that the bead location is within frame
    x_lb = min(X(:));
    x_ub = max(X(:));
    y_lb = min(Y(:));
    y_ub = max(Y(:));
    
    % set by size of bead and DNA strand
    z_lb = 0;
    z_ub = 600;
    
    
    % angles 
    theta_lb = 0;
    theta_ub = pi;
    phi_lb = -realmax;
    phi_ub = realmax;
    
    moment_lb = 0;
    moment_ub = 10*max(measured_fields(:))*(z_ub)^3;
    
    
    if options.Algorithm == 'levenberg-marquardt'
        lb = [];
        ub = [];
    else
        lb = [x_lb y_lb z_lb theta_lb phi_lb moment_lb] ./ rescale_vector;
        ub = [x_ub y_ub z_ub theta_ub phi_ub moment_ub] ./ rescale_vector;
    end
    
   
    
    
    
    [Tfit, ~, ~, ~, ~, ~, jacob] = lsqnonlin( @tmpFnt, T0, lb, ub, options);
    
    bead_parameter_output = Tfit .* rescale_vector;
    
    fisher_info = jacob'*jacob;
    
    if fixed_z_val
        xvals = [1 2 4 5 6];
        yvals = [1 2 4 5 6];
    
        fisher_info = fisher_info(xvals, yvals);
    end
    correlation_matrix = inv(fisher_info);
    
    
    fisher_info_precisions = (diag(correlation_matrix)).^(1/2);
    
    precisions = zeros(1,6);
    
    if fixed_z_val
        for i = [1,2,4,5,6]
            if i < 3
                precisions(i) = fisher_info_precisions(i)*rescale_vector(i);
            else
                precisions(i) = fisher_info_precisions(i-1)*rescale_vector(i);
            end
        end
        
    else
        
        for i = 1:6
            precisions(i) = fisher_info_precisions(i)*rescale_vector(i);
        end
    end





    function dy = tmpFnt(T)
        
        T(5) = mod(T(5), 2*pi);
        
%         theory = generateBeadFields(T.*rescale_vector, X, Y, nvDepth);
%         
%         temp = gradient_mask_simulated_b(theory, gradient_value, conversion_to_uT, step_size_nm);

        bead_params = T.*rescale_vector;
        
        if fixed_z_val
            bead_params(3) = initial_guess(3);
        end
        
        temp = simulate_bead_B_with_gradient_mask(bead_params, X, ...
            Y, nvDepth, pixel_side_length, subdivisions_per_axis, gradient_value, conversion_to_uT, ...
            front_binning_matrix, back_binning_matrix);
        
        
        
        if rescale_data_to_fit
           theory_mag = ( temp(:, :, 1).^2 + temp(:, :, 2).^2 + ...
               temp(:, :, 3).^2 ).^(1/2);
           data_mag = ( measured_fields(:, :, 1).^2 + measured_fields(:, :, 2).^2 + ...
               measured_fields(:, :, 3).^2 ).^(1/2);
           
           rescale_factor = max(data_mag(:)) / max(theory_mag(:));
           
           theory = temp*rescale_factor; 
           
        else 
            
            theory = temp;
        end
            
%% debug           
        if debug
            
            figure(10101);
            
            clf
            
            subplot(2,4,1);
            imagesc(squeeze(theory(:, :, 1)));
            pbaspect([1 1 1])
            colormap(linspecer)
            colorbar
            
            subplot(2,4,2);
            imagesc(squeeze(theory(:, :, 2)));
            pbaspect([1 1 1])
            colormap(linspecer)
            colorbar
            
            subplot(2,4,5);
            imagesc(squeeze(theory(:, :, 3)));
            pbaspect([1 1 1])
            colormap(linspecer)
            colorbar
            
            subplot(2,4,3);
            imagesc(squeeze(measured_fields(:, :, 1)));
            pbaspect([1 1 1])
            colormap(linspecer)
            colorbar
            
            subplot(2,4,4);
            imagesc(squeeze(measured_fields(:, :, 2)));
            pbaspect([1 1 1])
            colormap(linspecer)
            colorbar
            
            subplot(2,4,7);
            imagesc(squeeze(measured_fields(:, :, 3)));
            pbaspect([1 1 1])
            colormap(linspecer)
            colorbar
            
            drawnow
            
        end
        
%% stuff for masking

        mask_size = size(measured_fields);
    
        xx = 1:mask_size(1);
        yy = 1:mask_size(2);

        xx = xx - mean(xx(:));
        yy = yy - mean(yy(:));

        [mask_xgrid, mask_ygrid] = meshgrid(xx, yy);

        mask_rgrid = (mask_xgrid.^2 + mask_ygrid.^2).^(1/2);

        mask = (mask_rgrid < mask_max_dist) & (mask_rgrid > mask_min_dist);
%% generate differences            
            
        difference_x = squeeze(measured_fields(:, :, 1) - theory(:, :, 1));
        difference_y = squeeze(measured_fields(:, :, 2) - theory(:, :, 2));
        difference_z = squeeze(measured_fields(:, :, 3) - theory(:, :, 3));
        
        difference_x = difference_x./ squeeze(measured_field_precisions(:,:,1));
        difference_y = difference_y./ squeeze(measured_field_precisions(:,:,2));
        difference_z = difference_z./ squeeze(measured_field_precisions(:,:,3));
        
        
        masked_x = difference_x(mask);
        masked_y = difference_y(mask);
        masked_z = difference_z(mask);
        
        x_vals_for_lsq = masked_x;
        y_vals_for_lsq = masked_y; 
        z_vals_for_lsq = masked_z;
        
        %dy = [difference_x(mask), difference_y(mask), difference_z(mask)];
        
        
%        dy = x_vals_to_add(:) + y_vals_to_add(:) + z_vals_to_add(:);
        dy = x_vals_for_lsq(:) + y_vals_for_lsq(:);
%        dy = y_vals_to_add(:);
    end

end