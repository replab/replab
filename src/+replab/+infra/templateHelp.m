function res = templateHelp(templateName, element, helpCommand, strongIds, blackIds)
% Formats a documentation string for console output
%
% Args:
%   template (charstring): Name of the template to use, without path elements and without extension
%   element (`replab.infra.Element`): Element to expand the template with
%   helpCommand (charstring): Invocation of the help command, possibly including flags, without trailing space
%                             Examples would be 'help -f' or 'help'
%   strongIds (row cell vector of charstring): List of full identifiers to format with <strong> when supported)
%   plainIds (row cell vector of charstring): List of full identifiers to avoid linking
%
% Returns:
%   charstring: The interpreted documentation string
    codeBase = element.codeBase;
    templateFN = fullfile(codeBase.rootFolder, '+replab', '+infra', [templateName '.liquid']);
    template = replab.lobster.Template.load(templateFN);
    res = replab.infra.formatHelp(template.render(struct('el', element)), element, helpCommand, strongIds, blackIds);
end
