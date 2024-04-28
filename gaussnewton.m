function [x, y, a] = gaussnewton(omega, eta_values, Ic_list, alpha_list, x_tilde, y_tilde, a_tilde)
    max_iter = 100;
    tol = 1e-10;
    J = zeros(length(alpha_list), 3);
    new = [x_tilde; y_tilde; a_tilde];

    for iter = 1:max_iter
        old = new;
        
        F = (new(1) .* cos(new(2) .* cos(alpha_list) + new(3) .* sin(alpha_list))) - Ic_list;
        
        J(:, 1) = cos(new(2) .* cos(alpha_list) + new(3) .* sin(alpha_list));
        J(:, 2) = -new(1) .* cos(alpha_list) .* sin(new(2) .* cos(alpha_list) + new(3) .* sin(alpha_list));
        J(:, 3) = -new(1) .* sin(alpha_list) .* sin(new(2) .* cos(alpha_list) + new(3) .* sin(alpha_list));

        d = J \ (-F);
        new = old + d;

        if norm(old - new) < tol
            break;
        end
    end

    a = new(1) / eta_values;
    x = new(2) / omega;
    y = new(3) / omega;
end
