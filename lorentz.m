function f = lorentz(x0, gamma, I)
% This function returns a function handle to a Lorentzian function with the
% given parameters.

f =@handle;

function y = handle(x)
    y = ( (x-x0).^2 + gamma.^2 );
    y = I*gamma^2./y;
end
end