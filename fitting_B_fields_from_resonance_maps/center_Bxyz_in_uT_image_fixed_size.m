function [centered_Bxyz, sidelength, xrange, yrange] = ...
    center_Bxyz_in_uT_image_fixed_size(Bxyz_in_uT, input_sidelength, ...
    edge_cropping)

    debug = false;
    
    % for larger beads
    signal_threshold = 2.5; %uT
    
    % for turbobeads
    signal_threshold = 0.4; % uT
    
    
    Bx = Bxyz_in_uT(1+edge_cropping:end-edge_cropping,1+edge_cropping:end-edge_cropping,1);
    By = Bxyz_in_uT(1+edge_cropping:end-edge_cropping,1+edge_cropping:end-edge_cropping,2);
    Bz = Bxyz_in_uT(1+edge_cropping:end-edge_cropping,1+edge_cropping:end-edge_cropping,3);

%     Bx = Bxyz_in_uT(:,:,1);
%     By = Bxyz_in_uT(:,:,2);
%     Bz = Bxyz_in_uT(:,:,3);
    
    medfilt_size = 6;
    filtBx = medfilt2(Bx, [medfilt_size,medfilt_size]);
    filtBy = medfilt2(By, [medfilt_size,medfilt_size]);
    
    filt_B_mag = (filtBx.^2 + filtBy.^2).^(1/2);
    
    x = 1:size(filt_B_mag, 1);
    y = 1:size(filt_B_mag, 2);
    
    [X, Y] = meshgrid(x,y);
    
    CoM_normalization = sum(filt_B_mag(:));
    
    weighted_x = X.*filt_B_mag / CoM_normalization;
    weighted_y = Y.*filt_B_mag / CoM_normalization;
    
    xCoM = round(sum(weighted_x(:)));
    yCoM = round(sum(weighted_y(:)));
    
    is_above_threshold = filt_B_mag > signal_threshold;
    
    locations_where_above = find(is_above_threshold);

    if debug
        
        figure(19820); imagesc(is_above_threshold); pbaspect([ 1 1 1]);
        title('showing region centering on')
     
    end
    
    max_row = 1;
    min_row = size(is_above_threshold,2);
    
    max_col = 1;
    min_col = size(is_above_threshold, 1);
   
    
    for i = 1:length(locations_where_above)
        [row, col] = ...
            ind2sub(size(is_above_threshold), locations_where_above(i));
        
        if row < min_row
            min_row = row;
        end
        
        if row > max_row
            max_row = row;
        end
        
        if col < min_col
            min_col = col;
        end
        
        if col > max_col
            max_col = col;
        end
        
    end
    
    [topleft_col, topleft_row, sidelength] = ...
        generate_bounding_box(min_row, min_col, max_row, max_col, xCoM, yCoM);
    
    if topleft_col < 1
        topleft_col = 1;
    end
    
    if topleft_row < 1
        topleft_row = 1;
    end
    
    min_y = topleft_col;
    max_y = topleft_col + sidelength-1;
    
    min_x = topleft_row;
    max_x = topleft_row + sidelength-1;
    
    if sidelength < input_sidelength
        
        dif_in_sl = input_sidelength - sidelength;
        
        half_dif = floor(dif_in_sl/2);
        
        if mod(dif_in_sl, 2) == 0
            
            min_x = min_x - half_dif;
            min_y = min_y - half_dif;
            max_x = max_x + half_dif;
            max_y = max_y + half_dif;
            
        else
            
            min_x = min_x - half_dif;
            min_y = min_y - half_dif;
            max_x = max_x + half_dif+1;
            max_y = max_y + half_dif+1;
            
        end
        
    elseif sidelength > input_sidelength
        
        dif_in_sl = sidelength - input_sidelength ;
        
        half_dif = floor(dif_in_sl/2);
        
        if mod(dif_in_sl, 2) == 0
            
            min_x = min_x + half_dif;
            min_y = min_y + half_dif;
            max_x = max_x - half_dif;
            max_y = max_y - half_dif;
            
        else
            
            min_x = min_x + half_dif;
            min_y = min_y + half_dif;
            max_x = max_x - half_dif-1;
            max_y = max_y - half_dif-1;
            
        end
        
    end
            
    
    if min_x < 1 || min_y < 1 || max_x >  size(Bx, 1) || max_y >  size(Bx, 2)
        min_x = 1;
        min_y = 1;
        max_x = size(Bx, 1);
        max_y = size(Bx, 2);
        
        sidelength = size(Bx, 2);  
    end
        
    
    if sidelength > min(size(Bx, 1), size(Bx, 2))
        min_x = 1;
        min_y = 1;
        max_x = size(Bx, 1);
        max_y = size(Bx, 2);
        
        sidelength = size(Bx, 2);        
    end
    
    xrange = min_x:max_x;
    yrange = min_y:max_y;
    
    
    
    centered_Bxyz = zeros(length(xrange), length(yrange), 3);
    
    centered_Bxyz(:, :, 1) = Bx(xrange, yrange);
    centered_Bxyz(:, :, 2) = By(xrange, yrange);
    centered_Bxyz(:, :, 3) = Bz(xrange, yrange);
    
    
    