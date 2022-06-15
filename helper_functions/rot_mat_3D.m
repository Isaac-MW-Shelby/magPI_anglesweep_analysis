function [rotation_matrix] = rot_mat_3D(rotation_axis, rotation_angle)

    u = rotation_axis / norm(rotation_axis);
    
    ux = u(1);
    uy = u(2);
    uz = u(3);
    
    rotation_matrix = zeros(3,3);
    
    cos_t = cos(rotation_angle);
    sin_t = sin(rotation_angle);
    
    
    rotation_matrix(1,1) = cos_t + ux^2*(1-cos_t);
    rotation_matrix(1,2) = ux*uy*(1 - cos_t) - uz*sin_t;
    rotation_matrix(1,3) = ux*uz*(1 - cos_t) + uy*sin_t;
    
    rotation_matrix(2,1) = uy*ux*(1 - cos_t) + uz*sin_t;
    rotation_matrix(2,2) = cos_t + uy^2*(1-cos_t);
    rotation_matrix(2,3) = uy*uz*(1 - cos_t) - ux*sin_t;
    
    rotation_matrix(3,1) = uz*ux*(1 - cos_t) - uy*sin_t;
    rotation_matrix(3,2) = uz*uy*(1 - cos_t) + ux*sin_t;
    rotation_matrix(3,3) = cos_t + uz^2*(1-cos_t);
    
end