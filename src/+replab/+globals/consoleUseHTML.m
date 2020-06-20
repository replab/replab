function value = consoleUseHTML(newValue)
% Returns whether the console output accepts a subset of the HTML language
%
% The subset is the one parsed by MATLAB when using the GUI; it mostly accepts
% ``strong`` (for bold text) and ``a href=`` (for links) HTML tags.
    persistent storedValue
    if nargin == 1
        storedValue = newValue;
    end
    if isempty(storedValue)
        storedValue = ~replab.compat.isOctave && usejava('desktop');
    end
    value = storedValue;
end
