function [shifts] = spectrum_from_b_field(bfield,orientations)

    shifts = zeros(4,1);
    
    for i = 1:4
       shifts(i) = 2*dot(bfield, orientations(i, :)); 
    end
end

