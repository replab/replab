function str = lookupHelp(helpFunctionName, fullMode, identifier, variableName)
% Looks up the help contents for an object owned by RepLAB
%
% Args:
%   helpFunctionName (charstring): Name of the original help function
%   fullMode (logical): Whether to use the full display mode
%   identifier (charstring): Object identifier
%   variableName (charstring, optional): Variable name in the caller workspace
%
% Returns:
%   str: Help string
    parts = strsplit(identifier, '.');
    str = '';
    if length(parts) == 1 && nargin > 3
        str = [str sprintf('%s is an object of type %s.\n\n', firstPart, ...
                           replab.infra.repl.linkHelp(helpFunctionName, type, type))];
    end
    str = [str sprintf('--- help for %s ---\n\n', replab.infra.strong(identifier))];
    codeBase = replab.globals.codeBase;
    el = codeBase.get(parts{:});
    hasToggle = false;
    switch class(el)
      case 'replab.infra.Package'
        hasToggle = false;
        templateSuffix = 'package';
        docEl = el;
      case 'replab.infra.Class'
        hasToggle = true;
        templateSuffix = 'class';
        docEl = el;
      case 'replab.infra.Function'
        hasToggle = true;
        templateSuffix = 'function';
        docEl = el;
      case {'replab.infra.ConcreteClassElement', 'replab.infra.InheritedClassElement'}
        hasToggle = true;
        templateSuffix = 'class_element';
        docEl = el.declarations.findBestDocumented;
      otherwise
        error('replab:helpError', 'Object referenced by %s is of type %s', name, class(element));
    end
    if ~isempty(docEl) && (~isa(docEl, 'replab.infra.SourceElement') || docEl.doc.isempty)
        docEl = [];
    end
    if fullMode
        helpCurrent = [helpFunctionName ' -f'];
        helpToggled = helpFunctionName;
    else
        helpCurrent = helpFunctionName;
        helpToggled = [helpFunctionName ' -f'];
    end
    if hasToggle && fullMode
        templateName = ['help_full_' templateSuffix];
    else
        templateName = ['help_' templateSuffix];
    end
    fullId = el.fullIdentifier;
    str = [str replab.infra.repl.templateHelp(templateName, el, docEl, helpCurrent, {fullId}, {fullId})];
    if hasToggle && replab.globals.consoleUseHTML
        if fullMode
            link = replab.infra.repl.linkHelp(helpToggled, 'help mode', fullId);
            str = [str char(10) '   Toggle display to ' link char(10)];
        else
            link = replab.infra.repl.linkHelp(helpToggled, 'reference mode', fullId);
            str = [str char(10) '   Toggle display to ' link char(10)];
        end
    end
    if isa(el, 'replab.infra.SourceElement')
        sourceLink = replab.infra.repl.linkOpen('See source', '', el.absoluteFilename, el.startLineNumber);
        str = [str char(10) sourceLink char(10)];
    end
end
