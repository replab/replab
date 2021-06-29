function value = gapBinaryPath(newValue)
% Get/set the path to the GAP System binary
%
% (On a Linux system, this is the path to the ``gap.sh`` script).
%
% Args:
%   newValue (charstring or ``[]``, optional): Value to be stored if provided
%
% Returns:
%   integer: Stored integer value
    persistent storedValue % gets initialized to [] as per Matlab documentation
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        error('GAP system path not set');
    end
    value = storedValue;
end
