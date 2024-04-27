function x_opt = gyllensnitt(f,a,b,tol)
    % första itterationen
    update = @(a,b) (b-a)*(sqrt(5)-1)/2; % rest the update
    % update rule for x and y
    y = a + update(a,b);
    x = b - update(a,b);
    f1 = f(x);
    f2 = f(y);
    if f1 > f2
        a = x;
        prev = f1;
        prev_x = true;
    else
        b = y;
        prev = f2;
        prev_x = false;
    end

    while abs(b-a) > tol
       update = @(a,b) (b-a)*(sqrt(5)-1)/2; % rest update rule
       % update x and y
       if prev_x 
          x = y;
          y = a + update(a,b);
          f1 = prev;
          f2 = f(y);
       else
          y = x;
          x = b - update(a,b);
          f1 = f(x);
          f2 = prev;
       end
       % compare and cut down intervall
       if f1 > f2
           a = x;
           prev_x = true;
       else
           b = y;
           prev_x = false;
       end
    end
    x_opt = (a+b)/2;
end