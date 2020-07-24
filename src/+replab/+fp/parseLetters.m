function letters = parseLetters(word, names)
% Parses a word into letters
%
% Args:
%   word (charstring): Explicit word
%   names (cell(1,\*) of charstring): Generator names
%
% Returns:
%   integer(1,\*): Letters forming the word
    [ok, tokens] = replab.fp.Parser.lex(word, names);
    assert(ok, 'Unknown tokens in string');
    [pos, letters] = replab.fp.Parser.word(tokens, 1);
    assert(pos > 0, 'Malformed word');
    assert(tokens(1, pos) == replab.fp.Parser.types.END, 'Badly terminated word');
end
