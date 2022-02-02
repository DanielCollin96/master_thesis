function [param, residual, offset, slope] = fit_lorentz(data_x, data_y, x0, hwhm0, I0, offset0, slope0)
% Fit the data S by sum of Lorentzians using the positions of local maxima
% given in (x,y).  These values will be used for contructing reasonable
% starting values.
% Input:
%   data_x   --  x values of the data to be fit
%   data_y   --  y values of the data to be fit
%   x0       --  Initial values for the Lorentzian band position
%   hwhm0    --  Initial values for the HWHM
%   I0       --  Initial values for the intensities
%   offset0  --  Initial offset(s) for linear correction (vertical 
%            --  adjustments) (see below)
%   slope0   --  Initial slope for linear correction (see below)
%
% Output:
% CAUTION:  The behaviour of this function is governed by the number of
% requested output arguments:
%   [param,residual]
%          -- Plain Lorentzien Approximation
%   [param, residual, offset]
%          -- Lorentzian Approximation + vertical offset for each L.
%   [param, residual, offset, slope]
%          -- Lorentzian Approximation + global linear correction term
%
% The initial values for offset/slope are ignored if the corresponding
% output is not requested.

% NOTE:  We may experiment with setting bounds on the variables.  Sometimes
% we see negative values that result in useless and unstable Lorentzians
% NOTE:  We may want to clean up the Lorentzians in a run
% NOTE:  After we find no more residuals above the tolerance, we may do one
% final optimization for local improvement of the Lorentzians.  This could
% and should be combined with a clean up.

% Static enumeration of Lorentz-Approximation types.
ltypes.plain = 1;
ltypes.vertical = 2;
ltypes.linear = 3;

% Determine the type of approxmiation requested by caller.
if nargout == 2
    l_type = ltypes.plain;
    l_fun = @lorentz_plain;
elseif nargout == 3
    l_type = ltypes.vertical;
    l_fun = @lorentz_vertical;
elseif nargout == 4
    l_type = ltypes.linear;
    l_fun = @lorentz_linear;
else
    error('Wrong number of output arguments');
end

% The number of Lorentz function is encoded by the number of base points.
num_lo = length(x0);

% Default slope and offset
if nargin < 6 || isempty(offset0)
    if l_type == ltypes.vertical
        offset0 = zeros(num_lo, 1);
    else
        offset0 = 0;
    end
end

if nargin < 7 || isempty(slope0)
    slope0 = 0;
end

% Set the starting vector of the parameter vector for the Lorentzians
param = reshape([x0' ; hwhm0' ; I0'], 3*num_lo,1);

if l_type == ltypes.linear
    % The nonlinear variables vector is extended by two elements for the
    % linear correction.
    param = [param; offset0; slope0];
elseif l_type == ltypes.vertical
    % The nonlinear variables vector is extended by a vector of individual
    % offsets for each of the Lorentzian curves.
    param = [param; offset0];
end

% Set some options
options = optimset('Display', 'iter');
options.MaxFunEvals = 10000;
options.Jacobian = 'on';
options.TolFun = 1e-6;
options.MaxIter = 1000;

[param, ~, residual] = ...
    lsqcurvefit(l_fun, param, data_x, data_y, [], [], options);

% Reshape final parameters to more natural shape
if l_type == ltypes.linear
    offset = param(end-1);
    slope = param(end);
elseif l_type == ltypes.vertical
    offset = param(3*num_lo+1:end);
end

param = reshape(param(1:3*num_lo), 3, num_lo);
param = param';
end

function [y, J] = lorentz_plain(p,x)
y = evaluate_lorentz_sum(p, x);
if nargout > 1
    J = evaluate_lorentz_jacobian(p, x);
end

end % of function eval_lorentz_sum

function [y, J] = lorentz_vertical(p,x)
num_lo = idivide(length(p), int32(4));
p_lorentz = p(1:3*num_lo);
y = evaluate_lorentz_sum(p_lorentz, x);
y = y + p(3*num_lo+1:end);

if nargout > 1
    J = zeros(length(x), length(p));
    J(:,1:3*num_lo) = evaluate_lorentz_jacobian(p_lorentz, x);
    J(:,end) = 1;
end

end % of function eval_lorentz_sum


function [y, J] = lorentz_linear(p,x)
num_lo = idivide(length(p), int32(3));
p_lorentz = p(1:3*num_lo);
y = evaluate_lorentz_sum(p_lorentz, x);
y = y + p(end) * x + p(end-1);


if nargout > 1
    J = zeros(length(x), length(p));
    J(:,1:3*num_lo) = evaluate_lorentz_jacobian(p_lorentz, x);
    J(:,end-1) = 1;
    J(:,end) = x;
end

end % of function eval_lorentz_sum

function y = evaluate_lorentz_sum(p,x)
% Evaluate sum of Lorentzians.
num_lo = idivide(length(p), int32(3));
y = zeros(size(x));

for j=0:(num_lo - 1)
    x0 = p(j*3 + 1);
    gamma = p(j*3 + 2);
    I = p(j*3 + 3);
    
    y_tmp = (x - x0).^2 + gamma^2;
    
    y = y + I*gamma^2 ./ y_tmp;
    
end
end % of function eval_lorentz_sum

function J = evaluate_lorentz_jacobian(p,x)
% Evaluate Jacobian of a sum of Lorentzians

num_lo = idivide(length(p), int32(3));
J = zeros(length(x), length(p));
for j=0:(num_lo - 1)
    x0 = p(j*3 + 1);
    gamma = p(j*3 + 2);
    I = p(j*3 + 3);
    y_tmp = (x - x0).^2 + gamma^2;
    
    % partial derivative w.r.t. I
    J(:,j*3 + 3) = gamma^2 ./ y_tmp;
    
    % Update y_tmp to more convenient form
    y_tmp = y_tmp.^2;
    
    % partial derivative w.r.t. x0
    J(:, j*3 + 1) = 2*I*gamma^2 * (x - x0) ./ y_tmp;
    
    % partial derivative w.r.t. gamma
    J(:, j*3 + 2) = 2*gamma*I * (x - x0).^2 ./ y_tmp;
    
end

end % of function eval_lorentz_jacobian