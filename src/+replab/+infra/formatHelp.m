function res = formatHelp(txt, context, helpFunctionName, flags)
% Formats a documentation string for console output
%
% Args:
%   txt (charstring): String, possibly multiline (using ASCII 10 characters) 
%   context (`replab.infra.Element`): Element in the context of which the reference is interpreted
%
% Returns:
%   charstring: The interpreted documentation string
    lines = strsplit(txt, '\n');
    tab = sprintf('\t');
    tables = {}; % to memorize table lines
    tableStart = []; % toggles between [] (if not in a table), and the table start line
    for i = 1:length(lines)
        l = lines{i};
        
        %% Process Sphinx references
        ref_re = ['(?<!`)'       '`'        '(?!`)'];
        %        no ` before    backtick    no ` after
        parts = regexp(l, ref_re, 'split');
        txts = parts(1:2:end);
        refs = parts(2:2:end);
        for j = 1:length(refs)
            ref = refs{j};
            [el linkText] = replab.infra.resolve(ref, context, @(id) any(exist(id) == [3 4 5 6 8]));
            if isa(el, 'replab.Element')
                ref = replab.infra.linkHelp(helpFunctionName, linkText, el.fullIdentifier, flags);
            elseif isa(el, 'char')
                ref = replab.infra.linkHelp(helpFunctionName, linkText, el, flags);
            end
            refs{j} = ref;
        end
        
        parts(2:2:end) = refs;
        l = strjoin(parts);
        
        %% Replace Sphinx double backticks by single quotes
        l = strrep(l, '``', '''');        
        
        %% Identify tables
        if any(l == tab)
            if isempty(tableStart)
                tableStart = i;
            end
        else
            if ~isempty(tableStart)
                tables{end+1} = tableStart:(i-1);
                tableStart = [];
            end
        end
        
        %% Substitute formatted line
        lines{i} = l;
    end
    
    newLines = {};
    for i = length(tables):-1:1
        % process backwards so we don't crapify line numbers if the table formatting
        % ends up later taking a different number of lines
        table = tables{i};
        tableLines = lines(table);
    end
        
end
