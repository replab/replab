function value = defaultHelpFunction(newValue)
% Get/set the function handle corresponding to the default help function (one input/output variant)
%
% If called without argument, returns the stored value.
%
% The call with one argument changes the stored value; it should only be called
% as part of ``replab_init``.
%
% The stored value survives ``clear all`` as this function is locked.
%
% Args:
%   newValue (function_handle, optional): Sets the stored value
%
% Returns:
%   function_handle: Stored function handle
    mlock
    persistent storedValue
    if nargin == 1
        storedValue = newValue;
    end
    value = storedValue;
end
