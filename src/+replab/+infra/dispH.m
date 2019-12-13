function dispH(text, keyword, helpFunctionName, fullMode)
% Displays text while highlighting occurenced of keyword, and adding
% hyperlinks for words starting with "replab."

    if (replab.platformIsOctave) || (~usejava('desktop'))
        disp(text);
    else
        
        % We just do the processing by hand
        items = regexp(text, '(replab\.)[\w,\.]*', 'match');
        rest = regexp(text, '(replab\.)[\w,\.]*', 'split');
        for i = 1:length(items)
            matchedkeyword = regexp(items{i}, keyword);
            if ~isempty(matchedkeyword) && ((matchedkeyword(1) == 1) || isequal(items{i}(matchedkeyword(1)-1), '.') || isequal(items{i}(matchedkeyword(1)-1), ' '))
                %items{i} = regexprep(items{i}, keyword, ['<strong>', keyword, '</strong>']);
                items{i} = ['<strong>', items{i}, '</strong>'];
            else
                if fullMode
                    items{i} = ['<a href="matlab: ', helpFunctionName, '(''-f'',''', items{i}, ''')">', items{i}, '</a>'];
                else
                    items{i} = ['<a href="matlab: ', helpFunctionName, '(''', items{i}, ''')">', items{i}, '</a>'];
                end
            end
        end
        
        % We also boldify occurences of keyword which don't start with
        % replab.*
        for i = 1:length(rest)
            rest{i} = regexprep(rest{i}, [keyword, '(?!\w+)'], ['<strong>', keyword, '</strong>']);
        end
        
        % We recombine the text together
        text = '';
        for i = 1:length(items)
            text = [text, rest{i}, items{i}];
        end
        text = [text, rest{end}];
        
        disp(text);
    end
end
