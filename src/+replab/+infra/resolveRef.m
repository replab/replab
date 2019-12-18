function [el linkText] = resolveRef(context, ref, externalCheck)
% Resolves a Sphinx reference and returns the corresponding `.SourceElement`
%
% We accept two possible prefixes.
%
% The ``~`` prefix signifies that the link text is composed of the last element, for example
% ``~replab.Group.compose`` has link text ``compose`.
%   
% The ``.`` prefix signifies that we look up names first from the most precise scope then go
% to the general scope (thus from 4. to 1. below). 
%
% If the ``.`` prefix is not present, we look in the order 1. to 4.
%
% We search for the leading identifier in `ref` in the following scopes, from most general to most precise:
%
% 1. External objects: Objects outside the current `CodeBase`, for example in the Matlab global scope.
%    Those objects are looked up for using the `externalCheck` function handle.
%
% 2. Package names: the reference should start with a first-level package such as ``replab``.
%    We do not interpret second-level package names such as ``lobster`` (for ``replab.lobster``).
%
% 3. Class/function names in the current package
%
% 4. Method/property names in the current class (when applicable)
%
% When the leading identifier is of the first kind (1.), we do not resolve further, and simply return
% the charstring `ref` as the return value `el`.
%
% Otherwise, we return the corresponding `.SourceElement`.
%
% The returned link text will correspond to the reference, with the prefixes stripped; as an exception,
% when the ``~`` prefix is present, we return the rightmost identifier as the link text.
%
% Args:
%   context (`replab.infra.SourceElement`): Element in the context of which the reference is interpreted
%   ref (charstring): Reference
%   externalCheck (function_handle): Function that accepts a charstring identifier argument (without dots) and returns a logical value
%                                    that states whether the identifier is present in the external scope.
%  
% Returns
% -------
%   el: `.SourceElement` or ``[]``
%     The corresponding element, or ``[]`` if the reference could not be resolved
%   linkText: charstring
%     Extracted link text
    
end


% backticks for MATLAB declarations in other libraries -> help