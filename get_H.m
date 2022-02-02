function H = get_H(K,tspan,y0)
% This function solves the ODE y'=Ky, y(0)=y0, on the interval tspan, and
% returns the matrix H of relative concentrations.

odefun = @(t,y) K*y;
[~,H] = ode45(odefun, tspan, y0);
H = H';

end