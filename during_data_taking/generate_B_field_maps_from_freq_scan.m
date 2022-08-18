function [Bx, By, Bz] = generate_B_field_maps_from_freq_scan(resonance_map, input_signed_resonance_order)
%GENERATE_B_FIELD_MAPS_FROM_FREQ_SCAN Summary of this function goes here
%   Detailed explanation goes here


    squeezed_map = squeeze(resonance_map);
    
    sorted_map = sort(squeezed_map, 3);

    Z = sorted_map(:, :, 8:-1:5) - sorted_map(:, :, 1:4);
    
    resonance_order = abs(input_signed_resonance_order);
    resonance_signs = sign(input_signed_resonance_order);


    for res_num = 1:4
        
        res_ind = find(resonance_order == res_num);
        res_sign = resonance_signs(res_ind);
        
        Z(:, :, res_num) = squeeze(res_sign*squeezed_map(:, :, res_ind));
                
    end
    
    Bx = (sqrt(3))/8*(Z(:,:, 1) - Z(:, :, 2) - Z(:, :, 3) + Z(:, :, 4));
    By = (sqrt(3))/8*(Z(:,:, 1) - Z(:, :, 2) + Z(:, :, 3) - Z(:, :, 4));
    Bz = (sqrt(3))/8*(Z(:,:, 1) + Z(:, :, 2) - Z(:, :, 3) - Z(:, :, 4)); 
    
end

