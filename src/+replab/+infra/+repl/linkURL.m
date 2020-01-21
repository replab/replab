function str = linkURL(linkText, altText, url)
% Returns a HTML link that opens a source file in the editor if the console supports HTML, or alternative text
%
% Args:
%   linkText (charstring): Link text
%                         In there, ``%s`` is replaced by ``filename`` and ``%d`` by ``lineNumber``
%   altText (charstring): Alternative text if the console output does not support HTML
%                         In there, ``%s`` is replaced by ``filename`` and ``%d`` by ``lineNumber``
%   url (charstring): Linked URL
    if replab.settings.consoleUseHTML
        str = sprintf('<a href="matlab: web(''%s'', ''-browser'')">%s</a>', url, linkText);
    else
        str = altText;
    end
end
