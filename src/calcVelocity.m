function [v_vec, v_norm_vec] = calcVelocity(points, direction_vector, linear_velocity, rotation, angular_velocity, threshold)
    n_elements = size(points, 1);
    
    vcor = linear_velocity .* direction_vector .* 1000; 
    v_vec = ones(n_elements,1) .* vcor ;
    
    if rotation
        r_list = [points(:,1) points(:,2) points(:,3) + ones(n_elements,1) .* 100];
        v_sum = cross(ones(n_elements,1) .* angular_velocity, r_list) + v_vec;
        v_vec= round(v_sum, 15);
    end
    
    v_norm_vec = v_vec ./ vecnorm(v_vec, 2, 2);
    
    v_vec(abs(v_vec) < threshold) = 0;
    v_norm_vec(abs(v_norm_vec) < threshold) = 0;
end

