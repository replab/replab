function ref = processURLRef(ref)
% Process a Sphinx URL reference and returns the appropriate console text
%
% When using Matlab in the IDE, creates the relevant ``<a href='...'>`` reference, otherwise
% just prints the plain URL.
%
% Args:
%   ref (charstring): Sphinx URL reference of the form ``link text <URL ref>``.
%
% Returns:
%   charstring: Processed string
    tokens = regexp(ref, '^([^<]*)(<[^>]*>)$', 'tokens', 'once');
    if length(tokens) == 1
        tokens = {'' tokens{1}};
    end
    linkText = tokens{1};
    url = tokens{2};
    url = url(2:end-1); % strip the characters '<' and '>'
    ref = replab.infra.repl.linkURL(linkText, ref, url);
end
