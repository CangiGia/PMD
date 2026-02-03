function[value, isterminal, direction] = Ex_13_9_2(t,u)
    global R2

    value = R2 - norm(u(1:2) - u(4:5)); % gap between spheres
    isterminal = 1;
    direction = 1;
end
