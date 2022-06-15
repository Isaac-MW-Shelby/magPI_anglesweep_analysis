function [topleft_col, topleft_row, sidelength] = generate_bounding_box(min_row, min_col, max_row, max_col, xCoM, yCoM)



dist_com_to_max_col = max_col - xCoM;
dist_com_to_min_col = xCoM - min_col;

dist_com_to_max_row = max_row - yCoM;
dist_com_to_min_row = yCoM - min_row;

max_dist = max([dist_com_to_max_col, dist_com_to_min_col, ...
                dist_com_to_max_row, dist_com_to_min_row]);
            

sidelength = 2*max_dist;

topleft_col = xCoM - max_dist;
topleft_row = yCoM - max_dist;