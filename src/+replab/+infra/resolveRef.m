function [el linkText] = resolveRef(ref, context, isExternal)
% Resolves a Sphinx reference and returns the corresponding `.Element`
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
%    Those objects are looked up for using the `isExternal` function handle.
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
%   context (`replab.infra.Element`): Element in the context of which the reference is interpreted
%   ref (charstring): Reference
%   isExternal (function_handle): Function that accepts a charstring identifier argument (without dots) and returns a logical value
%                                 that states whether the identifier is present in the external scope.
%  
% Returns
% -------
%   el: `.Element` or ``charstring`` or ``[]``
%     The corresponding element, or ``[]`` if the reference could not be resolved
%   linkText: charstring
%     Extracted link text
    rx = ['^'  '([.~]*)' '([A-Za-z][A-za-z0-9_.]*)' '$'];
    %           prefix            identifier
    tokens = regexp(ref, rx, 'tokens', 'once');
    if length(tokens) == 1
        % For Octave compatibilty: octave will not separate an empty prefix
        tokens = {'' tokens{1}};
    end
    if length(tokens) ~= 2
        error('replab:resolveRefError', 'Invalid Sphinx reference %s', ref);
    end
    prefix = tokens{1};
    id = tokens{2};
    parts = strsplit(id, '.');
    head = parts{1};
    if isempty(parts) || any(cellfun(@isempty, parts))
        error('replab:resolveRefError', 'Invalid Sphinx reference %s', ref);
    end
    if ~ismember(prefix, {'', '.', '~', '~.', '.~'})
        error('replab:resolveRefError', 'Invalid prefix % in reference %s', prefix, ref);
    end
    if any(prefix == '~')
        linkText = parts{end};
    else
        linkText = id;
    end       
    if any(prefix == '.');
        order = [4 3 2 1];
    else
        order = [1 2 3 4];
    end
    headElement = [];
    for o = order
        switch o
          case 1 % lookup external scope
            if nargin >= 3 && ~isempty(isExternal) && isExternal(parts{1})
                el = id;
                return
            end
          case 2 % lookup absolute paths in code base
            e = context.codeBase.package(head);
            if ~isempty(e)
                headElement = e;
                break
            end
          case 3 % lookup current package
            if isa(context, 'replab.infra.Package')
                pkg = context;
            elseif isa(context, 'replab.infra.SourceElement')
                pkg = context.package;
            elseif isa(context, 'replab.infra.ClassElement')
                pkg = context.parentClass.package;
            else
                error('replab:resolveRefError', 'Invalid Sphinx reference %s', ref);
            end
            e = pkg.lookup(head);
            if ~isempty(e) && ~isa(e, 'replab.infra.Package')
                headElement = e;
                break
            end
          case 4 % lookup in current class if the context is a class
            if isa(context, 'replab.infra.Class')
                e = context.lookup(head);
                if ~isempty(e)
                    headElement = e;
                    break
                end
            end
        end
    end
    if isempty(headElement)
        el = [];
        return
    end
    
    % Circumventing direct assignment for octave
    if length(parts) == 1
        newParts = {};
    else
        newParts = parts(2:end);
    end
    el = headElement.get(newParts{:});
end
