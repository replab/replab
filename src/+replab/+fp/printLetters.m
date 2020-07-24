function word = printLetters(letters, names, times)
% Prints a word from its letters
%
% Args:
%   letters (integer(1,\*)): Letters of the word as generator indices
%   names (cell(1,\*) of charstring): Generator names
%   times (charstring, optional): Composition operator, default value ``' '``
    if nargin < 3
        times = ' ';
    end
    if isempty(letters)
        word = '1';
        return
    end
    word = '';
    sep = '';
    i = 1;
    letters = replab.fp.reduceLetters(letters);
    while i <= length(letters)
        e = sign(letters(i)); % exponent
        l = abs(letters(i));
        assert(e ~= 0);
        i = i + 1;
        while i <= length(letters) && abs(letters(i)) == l
            e = e + sign(letters(i));
            i = i + 1;
        end
        if e ~= 0
            word = [word sep names{l}];
            sep = times;
            if e ~= 1
                word = [word '^' num2str(e)];
            end
        end
    end
end
