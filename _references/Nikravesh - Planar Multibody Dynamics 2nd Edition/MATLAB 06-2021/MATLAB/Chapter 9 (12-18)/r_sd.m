function T_sd = r_sd(theta, theta_d, k, theta0, dc)
% This function computes the torque of a rotational spring-damper 

    T_sd = k*(theta - theta0) + dc*theta_d;
end
