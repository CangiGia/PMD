% Newton-Raphson iterative process
    x = input(' estimate for x = ?  ')
    results = zeros(15,3);
    
for n = 1:15
        f_x = x^3 -4*x^2 -44*x +96;
        df_dx = 3*x^2 - 8*x -44;
        results(n,:) = [n x f_x ];
    if abs(f_x) < 0.0001                % Is the constraint violated?
        break
    end
        d_x = -f_x/df_dx;               % Solve for the correction
        x = x + d_x;                    % Correct the estimate
end



