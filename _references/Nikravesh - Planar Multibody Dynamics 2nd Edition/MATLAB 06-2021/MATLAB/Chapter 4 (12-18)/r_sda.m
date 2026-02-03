function T_sda = r_sda(theta, theta_d, k, theta0, dc, Ta)
% This function computes the torque of a rotational combined 
%    spring-damper-actuator element 

    T_sda = k*(theta - theta0) + dc*theta_d + Ta;
end
