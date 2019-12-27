function tf = isOctave
% Returns true if and only if the platform is GNU/Octave.
%
% Taken from MOxUnit
%
% Returns:
%   logical: true if the platform on which this function is run is GNU/Octave, 
%            false otherwise (typically, if the platform is Matlab)
%
    persistent cached_tf

    if isempty(cached_tf)
        cached_tf = logical(exist('OCTAVE_VERSION', 'builtin'));
    end
    
    tf = cached_tf;
end
