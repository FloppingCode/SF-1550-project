function [x0_improved, y0_improved, a0_improved] = generate_guesses(Bound, omega, eta)
    g = Bound.un;

    ic1_diff = central_diff(Bound.x, Bound.y, g, pi/2, omega, Bound.s, 1e-8);
    ic2_diff = central_diff(Bound.x, Bound.y, g, 0, omega, Bound.s, 1e-8);
    is1 = simpson2(Bound.x, Bound.y, g, pi/2, omega, Bound.s);
    is2 = simpson2(Bound.x, Bound.y, g, 0, omega, Bound.s);
    ic1 = simpson(Bound.x, Bound.y, g, 0, omega, Bound.s);
    
    x0_improved = ic1_diff / (omega*is1);
    y0_improved =  -1 * (ic2_diff / (omega*is2));
    a0_improved = sqrt(ic1^2+is2^2) * 1/eta;
end

