function [p1, p2, p3, p4] = scatter4p3(Q, varargin)
% SCATTER4P3  Display quaternion array as four 3D scatter plots using
%             orthographic projection into four perpendicular 3-spaces.
%
% Takes the same parameters as the similarly-named Matlab functions, except
% that the first three parameters (X, Y, Z) are replaced by a single
% quaternion parameter. The varargin parameters are passed to the Matlab
% scatter3 function which is used to plot the four 3-dimensional plots, and
% must therefore conform to the requirements of that function.
%
% The output parameters are the plot axes to the four subplots, permitting
% the user to apply the set command to make adjustments to the axis limits
% and labels etc.

% Copyright (C) 2016, 2017 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% TODO The command axis equal affects only the last of the four subplots.
% Is there a way to make it operate on all four? Should the default be axis
% equal? Is there a way to set some status before plotting so that all axes
% have the same scale?

nargoutchk(0, 4)

for k = 1:length(varargin)
    if isa(varargin{k}, 'quaternion')
        error('Only the first parameter is permitted to be a quaternion.')
    end
end

if ~isvector(Q)
    error('Quaternion array must be a vector.')
    % TODO Why do we impose this limitation? We could access the parameter
    % using linear addressing, and thus cope with a matrix. NB The MATLAB
    % scatter functions insist on vectors, so at least we are consistent.
end

if ispure(Q)
    error(['Quaternion array must be a full quaternion: ',...
           'use scatter3 for pure quaternions'])
end

if ishold
    error('Hold is on, scatter4p3 does not currently support hold.')
    % TODO Consider whether this restriction could be removed, by turning
    % hold off, plotting the new data to each subplot, and ending with hold
    % on. This appears at first sight to be a complex problem to handle,
    % but it may be that a simple solution is possible with more thought.
end

% Find the extreme values in each of the four directions of 4-space, to
% determine the limits of the axes.

wmin = min(Q.w); wmax = max(Q.w);
xmin = min(Q.x); xmax = max(Q.x);
ymin = min(Q.y); ymax = max(Q.y);
zmin = min(Q.z); zmax = max(Q.z);

% Now derive vertical and horizontal maxima and minima so that the same
% scales are used vertically and horizontally in all four subplots.

vmin = min(xmin, zmin);
vmax = max(xmax, zmax);
hmin = min(wmin, ymin);
hmax = max(wmax, ymax);

% For the moment the vertical and horizontal scales are plotted equal, so
% find the minima and maxima of the four previous values.

plim = [min(vmin, hmin), max(vmax, hmax)];

% Four sub-plots are drawn, showing the scatter of values projected into
% four 3-spaces. Each 3-space omits one of the four standard axes (that's
% all we need to do to project the quaternion values into one of the
% standard 3-spaces).

p1 = subplot(2,2,1);
p2 = subplot(2,2,2);
p3 = subplot(2,2,3);
p4 = subplot(2,2,4);

scatter3(p1, Q.w, Q.x, Q.y, varargin{:}); % Omit Z.
scatter3(p2, Q.w, Q.x, Q.z, varargin{:}); % Omit Y.
scatter3(p3, Q.w, Q.y, Q.z, varargin{:}); % Omit X.
scatter3(p4, Q.x, Q.y, Q.z, varargin{:}); % Omit W.

set(p1, 'PlotBoxAspectRatio', [1 1 1]); % This ensures that the vertical
set(p2, 'PlotBoxAspectRatio', [1 1 1]); % and horizontal scales are equal.
set(p3, 'PlotBoxAspectRatio', [1 1 1]);
set(p4, 'PlotBoxAspectRatio', [1 1 1]);

set(p1, 'XLimMode', 'manual', 'YLimMode', 'manual', 'XLim', plim, 'Ylim', plim, 'Zlim', plim);
set(p2, 'XLimMode', 'manual', 'YLimMode', 'manual', 'Xlim', plim, 'YLim', plim, 'Zlim', plim);
set(p3, 'XLimMode', 'manual', 'YLimMode', 'manual', 'Xlim', plim, 'YLim', plim, 'Zlim', plim);
set(p4, 'XLimMode', 'manual', 'YLimMode', 'manual', 'Xlim', plim, 'YLim', plim, 'Zlim', plim);

set(get(p1, 'XLabel'), 'String', 'W'); % Label the axes, see above for the
set(get(p1, 'YLabel'), 'String', 'X'); % axis that is omitted in each case.
set(get(p1, 'ZLabel'), 'String', 'Y');

set(get(p2, 'XLabel'), 'String', 'W');
set(get(p2, 'YLabel'), 'String', 'X');
set(get(p2, 'ZLabel'), 'String', 'Z');

set(get(p3, 'XLabel'), 'String', 'W');
set(get(p3, 'YLabel'), 'String', 'Y');
set(get(p3, 'ZLabel'), 'String', 'Z');

set(get(p4, 'XLabel'), 'String', 'X');
set(get(p4, 'YLabel'), 'String', 'Y');
set(get(p4, 'ZLabel'), 'String', 'Z');

end

% $Id: scatter4p3.m 1004 2017-11-15 17:14:09Z sangwine $
