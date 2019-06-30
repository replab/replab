function test_version
% Test code to verify the Matlab or GNU Octave version.

% Copyright ? 2008, 2010, 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

% Notes on versions needed.
%
% From release 2 of QTFM versions of Matlab earlier than R2012a 7.14 are
% not accepted.
%
% The sylvester function, implemented from QTFM 2.5 requires Matlab version
% R2014a, since the Matlab function sylvester was not implemented in
% earlier versions. However, this requirement is superceded by the
% requirement for R2015b below.
%
% From release 2.4 of QTFM we require Matlab R2015b (8.6) or later, as we
% now implement the function numArgumentsFromSubscript, which resolves a
% long-standing issue with numel.


M = 8; % These are the minimum supported major and minor Matlab
N = 6; % version numbers.

disp('Checking Matlab version ...');

V = ver('Matlab'); % This gets a structure containing information about the
                   % Matlab version, if we are running under Matlab. If we
                   % are not, the result will be empty.
if ~isempty(V)

    S = V.Version; % This gets and displays a string containing the Matlab
    R = V.Release; % version and release. S starts with 7.7, 8.1 or similar.
                   % V will be a string such as '(R2015a)'.
                   
    disp([S, ' ', R]);       
    
    E = ['Matlab version must be ' num2str(M) '.' num2str(N) ' or later'];
    
    Dots = strfind(S, '.'); % Find the positions of the dots in the string.
    
    Major = str2double(S(1:Dots(1)-1));
    if length(Dots) == 1
        Minor = str2double(S(Dots(1) + 1:end));
    else
        Minor = str2double(S(Dots(1) + 1:Dots(2) - 1));
        % We ignore the sub-minor version if present. (e.g. 1 in 7.10.1).
    end
    
    if Major < 7
        error('Versions of Matlab earlier than 7 are not supported.')
    end
    
    if Major == M && Minor < N
        error(['Matlab version must be ' ...
               num2str(M) '.' num2str(N) ' or later.']);
    end
    
    % Otherwise Major must be greater than M and we are OK, without
    % checking the minor version.
    
    disp('Matlab version is OK.');

else
    
    disp('Not running under Matlab, checking for Octave ...');
    
    V = ver('Octave');
    if isempty(V)
        error('Not running under Matlab or Octave, giving up!');
    end
    
    S = V.Version;
    if S(2) == '.' && S(4) == '.'
        if str2double(S(1)) < 3
           error('Minimum major version required is 3.');
        else
            if str2double(S(1)) == 3 && str2double(S(3)) < 2
                error('Minimum minor version required is 2.');
            end
        end
    else
        error('Cannot make sense of Octave version string.')
    end

    disp('Octave version is OK.');
    
end

% $Id: test_version.m 1004 2017-11-15 17:14:09Z sangwine $
