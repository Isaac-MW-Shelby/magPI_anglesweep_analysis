function [centered_Bxyz, sidelength, xrange, yrange] = ...
    center_Bxyz_in_uT_image_fixed_size(Bxyz_in_uT, input_sidelength, ...
    edge_cropping)

    debug = true;
    
    % for larger beads
    signal_threshold = 2.5; %uT
    
    % for turbobeads
    signal_threshold = 0.6; % uT
    
    
    
    Bx = Bxyz_in_uT(1+edge_cropping:end-edge_cropping,1+edge_cropping:end-edge_cropping,1);
    By = Bxyz_in_uT(1+edge_cropping:end-edge_cropping,1+edge_cropping:end-edge_cropping,2);
    Bz = Bxyz_in_uT(1+edge_cropping:end-edge_cropping,1+edge_cropping:end-edge_cropping,3);

    xy_size = size(Bx,1);
    
    medfilt_size = 6;
    
    B_mag = (Bx.^2 + By.^2 + Bz.^2).^2;
    
    filt_B_mag = medfilt2(B_mag, [medfilt_size, medfilt_size]);
    
    x = 1:size(filt_B_mag, 1);
    y = 1:size(filt_B_mag, 2);
    
    [X, Y] = meshgrid(x,y);
    
    is_above_threshold = filt_B_mag > signal_threshold;

    if debug
        
        figure(19820); imagesc(is_above_threshold); pbaspect([ 1 1 1]);
        title('showing region centering on')
     
    end
    
    number_non_zero = sum(is_above_threshold(:));
    
    weighted_x = X.*is_above_threshold / number_non_zero;
    weighted_y = Y.*is_above_threshold / number_non_zero;
    
    xCoM = round(sum(weighted_x(:)));
    yCoM = round(sum(weighted_y(:)));
   
    
    min_x = xCoM - floor(input_sidelength/2);
    max_x = xCoM + ceil(input_sidelength/2);
    
    min_y = yCoM - floor(input_sidelength/2);
    max_y = yCoM + ceil(input_sidelength/2);
    
    if input_sidelength >= xy_size
        
        min_x = 1;
        min_y = 1;
        
        max_x = xy_size;
        max_y = xy_size;
        
    end
    
    
    if min_x < 1
       
        difference = 1 - min_x;
        
        min_x = min_x + difference;
        max_x = max_x + difference;
        
    end
    
    if min_y < 1
       
        difference = 1 - min_y;
        
        min_y = min_y + difference;
        max_y = max_y + difference;
        
    end
    
    if max_x > xy_size
        
        difference = max_x - xy_size;
        
        max_x = max_x - difference;
        min_x = min_x - difference;
        
    end
    
    if max_y > xy_size
        
        difference = max_y - xy_size;
        
        max_y = max_y - difference;
        min_y = min_y - difference;
        
    end
    
    if min_y < 1
        min_y = 1;
    end
    
    if min_x < 1
        min_x = 1;
    end
    
%     xrange = min_x:max_x;
%     yrange = min_y:max_y;
    
    
    sidelength = input_sidelength;
    
    
    xrange = min_y:max_y;
    yrange = min_x:max_x;
    
    
    
    centered_Bxyz = zeros(length(xrange), length(yrange), 3);
    
    centered_Bxyz(:, :, 1) = Bx(xrange, yrange);
    centered_Bxyz(:, :, 2) = By(xrange, yrange);
    centered_Bxyz(:, :, 3) = Bz(xrange, yrange);
    
    
    