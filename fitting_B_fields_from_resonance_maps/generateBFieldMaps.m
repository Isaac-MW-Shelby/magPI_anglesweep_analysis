function [per_pixel_B_fields, per_pixel_B_precisions] = ...
    generateBFieldMaps(shift_data, resonance_assignment)


per_pixel_B_fields = zeros(size(shift_data,1), size(shift_data,2), size(shift_data,3), 3);
per_pixel_B_precisions = zeros(size(per_pixel_B_fields));


for i = 1:size(shift_data, 1)
    
    resonance_signs = sign(resonance_assignment(i,:));
    resonance_order = abs(resonance_assignment(i, :));
    
    Z = zeros(size(squeeze(shift_data(i, :, :, :))));
    
    for res_num = 1:4
        
        res_ind = find(resonance_order == res_num);
        res_sign = resonance_signs(res_ind);
        
        Z(:, :, res_num) = squeeze(res_sign*shift_data(i, :, :, res_ind));
                
    end
    
    Bx = (sqrt(3))/8*(Z(:,:, 1) - Z(:, :, 2) - Z(:, :, 3) + Z(:, :, 4));
    By = (sqrt(3))/8*(Z(:,:, 1) - Z(:, :, 2) + Z(:, :, 3) - Z(:, :, 4));
    Bz = (sqrt(3))/8*(Z(:,:, 1) + Z(:, :, 2) - Z(:, :, 3) - Z(:, :, 4)); 
    
    per_pixel_B_fields(i, :, :, 1) = Bx;
    per_pixel_B_fields(i, :, :, 2) = By;
    per_pixel_B_fields(i, :, :, 3) = Bz;
 
    
end

end


