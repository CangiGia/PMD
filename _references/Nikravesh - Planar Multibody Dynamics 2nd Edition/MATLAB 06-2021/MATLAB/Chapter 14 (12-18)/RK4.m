function [T, zT] = RK4(fun, Tspan, z)
% 4th order Runge-Kutta algorithm
    T = Tspan'; zT = z'; nT = size(Tspan);
for i=1:(nT(2) - 1)
    t  = Tspan(i); dt = Tspan(i + 1) - Tspan(i);
    f1 = fun(t, z);
    f2 = fun((t + dt/2), (z + dt*f1/2));
    f3 = fun((t + dt/2), (z + dt*f2/2));
    f4 = fun((t + dt  ), (z + dt*f3  ));
    z  = z + dt*(f1 + 2*(f2 + f3) + f4)/6;
    zT = [zT; z'];
end
end
