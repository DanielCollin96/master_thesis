function [p_mat, offset, slope] = greedy_lorentz(data_x, data_y, debug)
% A simple greedy heuristic for fitting data with Lorentzians.
%
% Paramters
%   data_x    x values of the data to be fit
%   data_y    y values of the data to be git
%   debug     flag for visual debugging
%
% CAUTION:  The behaviour of this function is governed by the number of
% output arguments requested.  
%   [p_mat]
%          -- Plain Lorentzien Approximation
%   [p_mat, offset]
%          -- Lorentzian Approximation + vertical offset for each L.
%   [p_mat, offset, slope]
%          -- Lorentzian Approximation + global linear correction term
%
% Output description:
%    p_mat    k x 3 matrix where every gives a triple (x0, gamma, I) of
%             parameters describing a Lorentzian.  'k' corresponds to the
%             number of Lorentzians.
%
%    offset   offset on the y-axis of the linear correction (or individual
%             offsets for each of the Lorentzians in p_mat.
%    slope    slope of the linear correction

% Static enumeration of Lorentz-Approximation types.
ltypes.plain = 1;
ltypes.vertical = 2;
ltypes.linear = 3;

if nargin < 3 || isempty(debug)
    debug = false;
end

if nargout == 1
    l_type = ltypes.plain;
elseif nargout == 2
    l_type = ltypes.vertical;
elseif nargout == 3
    l_type = ltypes.linear;
else
    error('Incompatible number of output arguments.');
end

eps_drop_I = 100;

[maxres_val, maxres_ind] = max(data_y);

x0 = data_x(maxres_ind);
I0 = maxres_val;
hwhm0 = 1;
offset0 = 1;
slope0 = 0;

while (maxres_val > eps_drop_I)
    if l_type == ltypes.plain
        [p_mat, residual] = fit_lorentz(data_x, data_y, x0, hwhm0, I0);
    elseif l_type == ltypes.linear
        [p_mat, residual, offset, slope] = ...
            fit_lorentz(data_x, data_y, x0, hwhm0, I0, offset0, slope0);
    elseif l_type == ltypes.vertical
        [p_mat, residual, offset] = ...
            fit_lorentz(data_x, data_y, x0, hwhm0, I0, offset0);
    else
        % Cannot happen
        error('fixme');
    end
    
    [maxres_val, maxres_ind] = min(residual);
    maxres_val = -maxres_val;

    if debug
        fprintf('Maximum residual (mean): %f (%f)\n', maxres_val,...
            mean(abs(residual)));
        if l_type == ltypes.linear
            plot_lorentz_fit(data_x, data_y, p_mat, 0, offset, slope);
        elseif l_type == ltypes.plain
            plot_lorentz_fit(data_x, data_y, p_mat);
        elseif l_type == ltypes.vertical
            plot_lorentz_fit(data_x, data_y, p_mat, 0, offset);
        else
            % Cannot happen
            error('fixme');
        end
        hold on
        plot(data_x(maxres_ind), maxres_val, 'ro', 'MarkerSize', 20,...
            'MarkerFaceColor', 'c');
        hold off
        pause;
    end
    x0 = [p_mat(:,1); data_x(maxres_ind)];
    hwhm0 = [p_mat(:,2); 1];
    I0 = [p_mat(:,3); maxres_val];    
    
    if l_type == ltypes.vertical
        offset0 = [offset; 1];
    elseif l_type == ltypes.linear
        offset0 = offset;
        slope0 = slope;
    elseif l_type == ltypes.plain
        % No initial values for slope or offset needed
    else
        % Cannot happen
        error('fixme');
    end
end

if debug
    hold off
end

end % of greedy_lorentz