function dispH(text, fullKeyword, helpFunctionName, fullMode)
% Display with highlight formatting
% 
% Displays text while highlighting occurenced of keyword, and adding
% hyperlinks for words starting with "replab." (which are different than
% the highlighted keyword).
%
% Currently, highlighting and hyperlinks are only available on matlab in
% desktop mode.
% 
% Args:
%     text (charstring): text to be displayed
%     fullKeyword (charstring): fully qualified class/function name to be
%         highlighted in the string
%     helpFunctionName (charstring): name of the "help" function
%     fullMode (bool): whether links should be created to the help in
%         simple or full mode

    if (replab.platformIsOctave) || (~usejava('desktop'))
        disp(text);
    else
        % extract the final element
        keyword = strsplit(fullKeyword, '.');
        keyword = keyword{end};
        
        % We just do the processing by hand
        items = regexp(text, '(replab\.)[\w,\.]*', 'match');
        rest = regexp(text, '(replab\.)[\w,\.]*', 'split');
        for i = 1:length(items)
            if isequal(items{i}, fullKeyword)
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
