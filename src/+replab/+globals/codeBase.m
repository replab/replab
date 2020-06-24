function value = codeBase(newValue)
% Get/set the parsed code used in the help system
%
% If called without argument, returns the stored value.
%
% The call with one argument changes the stored value; it should only be called
% as part of ``replab_init``.
%
% The stored value survives ``clear all`` as this function is locked.
%
% Call ``replab.globals.codeBase([])`` or ``help -clear`` to clear it.
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
        if isempty(newValue)
            value = [];
            return
        end
    end
    if isempty(storedValue)
        disp('Building code index');
        storedValue = replab.infra.crawl(fullfile(replab.globals.replabPath, 'src'));
    end
    value = storedValue;
end
