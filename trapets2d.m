function [integral] = trapets2d(func, x_min, x_max, y_min, y_max, num_points, a)
    dx = (x_max - x_min) / num_points;
    dy = (y_max - y_min) / num_points;
    integral = 0;  
    for i = 1:num_points
        for j = 1:num_points
            x = x_min + (i - 1) * dx;
            y = y_min + (j - 1) * dy;
            integral = integral + func(x, y, a); 
        end
    end
    integral = integral * dx * dy;
end
