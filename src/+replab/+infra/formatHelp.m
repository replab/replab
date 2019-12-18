function res = formatHelp(txt, context, helpFunctionName, flags)
% Formats a documentation string for console output
%
% Args:
%   txt (charstring): String, possibly multiline (using ASCII 10 characters) 
%   context (`replab.infra.Element`): Element in the context of which the reference is interpreted
%
% Returns:
%   charstring: The interpreted documentation string
    lines = strsplit(txt, '\n', 'CollapseDelimiters', false);
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
            [el linkText] = replab.infra.resolveRef(ref, context, @(id) any(exist(id) == [3 4 5 6 8]));
            if isa(el, 'replab.infra.Element')
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
    
    %% Format identified tables
    for i = length(tables):-1:1
        % process backwards so we don't crapify line numbers if the table formatting
        % ends up later taking a different number of lines
        tableLN = tables{i};
        nRows = length(tableLN);
        tableLines = lines(tableLN);
        nCols = max(cellfun(@(l) sum(l == tab), tableLines)) + 1;
        table = cell(nRows, nCols);
        for r = 1:nRows
            cols = strsplit(tableLines{r}, '\t');
            table(r, 1:length(cols)) = cols;
        end
        tableLines = replab.infra.align(table, repmat('l', 1, nCols));
        lines = {lines{1:tableLN(1)-1} tableLines{:} lines{tableLN(end)+1:end}};
    end
    
    res = strjoin(lines, '\n');
end
