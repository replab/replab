function res = templateHelp(templateName, el, docEl, helpCommand, strongIds, blackIds)
% Formats a documentation string for console output
%
% Args:
%   template (charstring): Name of the template to use, without path elements and without extension
%   el (`replab.infra.Element`): Element passed as the parameter ``el`` in the template context
%   docEl (`replab.infra.Element` or ``[]``): Element containing a documentation comment block. Can differ
%                                             from ``el`` for class methods/properties, when the documentation is
%                                             inherited. Note that the rendered template references are expanded with
%                                             respect to ``docEl``, not ``el``, when ``docEl`` is present; so care
%                                             must be taken in template output to write full identifiers in references.
%                                             If no documentation is available, then ``docEl`` is ``[]``.
%   helpCommand (charstring): Invocation of the help command, possibly including flags, without trailing space
%                             Examples would be 'help -f' or 'help'
%   strongIds (row cell vector of charstring): List of full identifiers to format with <strong> when supported)
%   plainIds (row cell vector of charstring): List of full identifiers to avoid linking
%
% Returns:
%   charstring: The interpreted documentation string
    codeBase = el.codeBase;
    templateFN = fullfile(codeBase.rootFolder, '+replab', '+infra', [templateName '.liquid']);
    template = replab.lobster.Template.load(templateFN);
    context = struct('el', el, 'docEl', docEl);
    if isempty(docEl)
        renderEl = el;
    else
        renderEl = docEl;
    end
    res = replab.infra.formatHelp(template.render(context), renderEl, helpCommand, strongIds, blackIds);
end
