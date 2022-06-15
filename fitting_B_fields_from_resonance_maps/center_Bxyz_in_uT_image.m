function [centered_Bxyz, sidelength] = center_Bxyz_in_uT_image(Bxyz_in_uT)
    

    signal_threshold = 2.5; %uT

    
    Bx = Bxyz_in_uT(:,:,1);
    By = Bxyz_in_uT(:,:,2);
    Bz = Bxyz_in_uT(:,:,3);
    
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
    
%     figure(19820); imagesc(is_above_threshold); pbaspect([ 1 1 1]);
%     title('showing region centering on')
    
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
    

%     figure(19821);
%     colormap(linspecer);
% 
%     subplot(2,3,1);
%     imagesc(Bx);
%     title("raw Bx");
%     axis square;
%     colorbar;
%     set(get(colorbar,'label'),'String','Bx (\muT)');
%     
%     subplot(2,3,2);
%     imagesc(By);
%     title("raw By");
%     axis square;
%     colorbar;
%     set(get(colorbar,'label'),'String','By (\muT)');
% 
%     subplot(2,3,3);
%     imagesc(Bz);
%     title("raw Bz");
%     axis square;
%     colorbar;
%     set(get(colorbar,'label'),'String','Bz (\muT)');
%     
%     subplot(2,3,4);
%     imagesc(centered_Bxyz(:,:,1));
%     title("centered Bx");
%     axis square;
%     colorbar;
%     set(get(colorbar,'label'),'String','Bx (\muT)');
%     
%     subplot(2,3,5);
%     imagesc(centered_Bxyz(:,:,2));
%     title("centered By");
%     axis square;
%     colorbar;
%     set(get(colorbar,'label'),'String','By (\muT)');
% 
%     subplot(2,3,6);
%     imagesc(centered_Bxyz(:,:,3));
%     title("centered Bz");
%     axis square;
%     colorbar;
%     set(get(colorbar,'label'),'String','Bz (\muT)');
    
    
end