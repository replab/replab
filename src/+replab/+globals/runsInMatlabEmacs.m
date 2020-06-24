function value = runsInMatlabEmacs(newValue)
% Returns whether the MATLAB/Octave shell is running in the matlab-emacs mode
    persistent storedValue
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = false;
    end
    value = storedValue;
end
