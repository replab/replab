function tf = platformIsOctave
% Returns true if and only if the platform is GNU/Octave.
%
% tf = replab.platformIsOctave()
%
% Output:
%   tf      true if the platform on which this function is run is
%           GNU/Octave, false otherwise (typically, if the platform is
%           Matlab)
%
% Taken from MOxUnit

    persistent cached_tf;
    
    if islogical(cached_tf)
        tf = cached_tf;
        return;
    end

    tf = logical(exist('OCTAVE_VERSION', 'builtin'));
    cached_tf = tf;
end


