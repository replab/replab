function letters = parseLetters(word, names)
    [ok, tokens] = replab.fp.Parser.lex(word, names);
    assert(ok, 'Unknown tokens in string');
    [pos, letters] = replab.fp.Parser.word(tokens, 1);
    assert(pos > 0, 'Malformed word');
    assert(tokens(1, pos) == replab.fp.Parser.types.END, 'Badly terminated word');
end
