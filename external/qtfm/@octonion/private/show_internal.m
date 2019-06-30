function show_internal(name, value)
% Displays the eight components of an octonion (array). This code is called
% from show.m and displayall.m, which must pass the name of the variable as
% a string. (It may be empty if the value is anonymous, as result for
% example, from the call 'show(qi)'.)

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

if isempty(value)
    
    if ~isempty(name)
        disp(' ');
        disp([name, ' =']);
    end

    % There are no numeric values to be output, so we simply output a
    % description of the empty value, depending on its size, which may be
    % 0-by-0, or 0-by-N, and it may be a matrix (2-dimensional) or an array
    % (more than two-dimensional). We can output something informative in
    % all these cases and not just '[]'.
    
    S = blanks(5);
    d = size(value);
    if sum(d) == 0
        S = [S '[] octonion'];
    else
        S = [S 'Empty octonion'];
        l = length(d);
        if l == 2
            S = [S ' matrix: '];
        else
            S = [S ' array: '];
        end
        for k = 1:l
            S = [S, num2str(d(k))];
            if k == l
                break % If we have just added the last dimension, no need
                      % for another multiplication symbol.
            end
            S = [S, '-by-'];
        end
    end
    disp(S)
else
    % The octonion is not empty, therefore we can output numeric data. We
    % do this by outputting the seven or eight components one by one. There
    % is a special case if the octonion is scalar, since we can output
    % this on one line, using disp.
    
    if isscalar(value)
        disp(' ');
        disp([name, ' =']);
        disp(' ');
        disp(value);
        disp(' ');
    else
        disp(' ');
        % The choice of labels here is tricky. In the quaternion case, the
        % labels are SXYZ, and it seems sensible to use the same for the
        % first four components of the octonion. What to use for the others
        % is difficult. ABCD in some sense 'follows on' from XYZ. In an
        % earlier version of this code, the labels used were STUVWXYZ,
        % which has some merit, but it meant that the 2nd to 4th components
        % were labelled differently to the quaternion case, and V was used
        % for one of them, when v() is used for the vector part.
        Labels = 'SXYZABCD';
        for n = 1:8
            if n == 1 && ispure(value)
              continue % Skip output of an empty scalar part.
            end
            L = Labels(n);
            if isempty(name)
                disp([L ' =']); disp(' ');
            else
                disp([name '.' L ' =']); disp(' ');
            end
            disp(component(value, n));
        end
    end
end

% $Id: show_internal.m 1004 2017-11-15 17:14:09Z sangwine $
