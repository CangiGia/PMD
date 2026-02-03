function T_s = r_s(theta, k, theta0)
% This function computes the torque of a rotational spring 

    T_s = k*(theta - theta0);
end
