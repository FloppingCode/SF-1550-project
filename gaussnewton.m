function [x0,y0,a, counter] = gaussnewton(f_c, Ic_list, aa_list, w, x_tilde, y_tilde, a_tilde)

    new = [w * x_tilde; w * y_tilde; f_c * a_tilde];
    iterations = 1e6;
    tolerance = 1e-10;
    counter = 0;

    for i = 1:iterations
        old = new;
        I_c = Func(old, aa_list, Ic_list);
        J = jacobian(old(1), old(2), old(3), aa_list);
        d = J \ (-I_c);

        new = old + d;
        counter = counter + 1;
        if norm(old - new) < tolerance
            break
        end
    end

    x0 = new(1, 1) / w;
    y0 = new(2, 1) / w;
    a = new(3, 1) / f_c;

end

function I_c = Func(old, aa_list, Ic_list)
    x_tilde = old(1);
    y_tilde = old(2);
    a_tilde = old(3);
    I_c = a_tilde .* cos(x_tilde .* cos(aa_list) + y_tilde .* sin(aa_list)) - Ic_list;
end

function J = jacobian(x_tilde, y_tilde, a_tilde, aa_list)
    n = length(aa_list);
    J = zeros(n, 3);
    for i = 1:n
        J(i, 1) = -a_tilde * sin(x_tilde * cos(aa_list(i)) + y_tilde * sin(aa_list(i))) * cos(aa_list(i));
        J(i, 2) = -a_tilde * sin(x_tilde * cos(aa_list(i)) + y_tilde * sin(aa_list(i))) * sin(aa_list(i));
        J(i, 3) = cos(x_tilde * cos(aa_list(i)) + y_tilde * sin(aa_list(i)));
    end
end
