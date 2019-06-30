function scatter4(Q, varargin)
% SCATTER4  Display quaternion array as 4D scatter plot using orthographic
%           projection onto six perpendicular planes.
%
% Takes the same parameters as the similarly-named Matlab functions, except
% that the first three parameters (X, Y, Z) are replaced by a single
% quaternion parameter. The varargin parameters are passed to the Matlab
% scatter function which is used to plot the six 2-dimensional plots, and
% must therefore conform to the requirements of that function.

% TODO Consider adding output parameters as per scatter4p3.

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

nargoutchk(0, 0)

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
    error('Hold is on, scatter4 does not currently support hold.')
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
% scales are used vertically and horizontally in all six subplots.

vmin = min(xmin, zmin);
vmax = max(xmax, zmax);
hmin = min(wmin, ymin);
hmax = max(wmax, ymax);

% For the moment the vertical and horizontal scales are plotted equal, so
% find the minima and maxima of the four previous values.

plim = [min(vmin, hmin), max(vmax, hmax)];

% Six sub-plots are drawn, showing the scatter of values projected onto six
% planes: (1, i), (1, j), (1, k) and (i, j), (j, k) and (i, k).

p1 = subplot(2,3,1); p2 = subplot(2,3,2);
p3 = subplot(2,3,3); p4 = subplot(2,3,4);
p5 = subplot(2,3,5); p6 = subplot(2,3,6);

scatter(p1, Q.w, Q.x, varargin{:});
scatter(p2, Q.w, Q.y, varargin{:});
scatter(p3, Q.w, Q.z, varargin{:});
scatter(p4, Q.x, Q.y, varargin{:});
scatter(p5, Q.y, Q.z, varargin{:});
scatter(p6, Q.z, Q.x, varargin{:});

set(p1, 'PlotBoxAspectRatio', [1 1 1]); % This ensures that the vertical
set(p2, 'PlotBoxAspectRatio', [1 1 1]); % and horizontal scales are equal.
set(p3, 'PlotBoxAspectRatio', [1 1 1]);
set(p4, 'PlotBoxAspectRatio', [1 1 1]);
set(p5, 'PlotBoxAspectRatio', [1 1 1]);
set(p6, 'PlotBoxAspectRatio', [1 1 1]);

set(p1, 'XLimMode', 'manual', 'YLimMode', 'manual', 'XLim', plim, 'Ylim', plim);
set(p2, 'XLimMode', 'manual', 'YLimMode', 'manual', 'Xlim', plim, 'YLim', plim);
set(p3, 'XLimMode', 'manual', 'YLimMode', 'manual', 'Xlim', plim, 'YLim', plim);
set(p4, 'XLimMode', 'manual', 'YLimMode', 'manual', 'Xlim', plim, 'YLim', plim);
set(p5, 'XLimMode', 'manual', 'YLimMode', 'manual', 'Xlim', plim, 'YLim', plim);
set(p6, 'XLimMode', 'manual', 'YLimMode', 'manual', 'Xlim', plim, 'YLim', plim);

set(get(p1, 'XLabel'), 'String', 'W');
set(get(p1, 'YLabel'), 'String', 'X');
set(get(p2, 'XLabel'), 'String', 'W');
set(get(p2, 'YLabel'), 'String', 'Y');
set(get(p3, 'XLabel'), 'String', 'W');
set(get(p3, 'YLabel'), 'String', 'Z');
set(get(p4, 'XLabel'), 'String', 'X');
set(get(p4, 'YLabel'), 'String', 'Y');
set(get(p5, 'XLabel'), 'String', 'Y');
set(get(p5, 'YLabel'), 'String', 'Z');
set(get(p6, 'XLabel'), 'String', 'Z');
set(get(p6, 'YLabel'), 'String', 'X');

end

% $Id: scatter4.m 1004 2017-11-15 17:14:09Z sangwine $
