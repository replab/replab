classdef Parser
% Parses a group presentation, or a group element expressed as a product of generators
%
% The grammar is given by:
%
%   <presentation> ::= '<' <generator> (','? <generator>)* '|' <relators> '>'
%   <relators> ::= (<equation> (',' <equation>)*)?
%   <equation> ::= <nonEmptyWord> ('=' <nonEmptyWord>)*
%   <word> ::= <nonEmptyWord>?
%   <nonEmptyWord> ::= <component> ('*' <component> | '/' <component> | <component>)*
%   <component> ::= <part> <exponent>?
%   <part> ::= <identity> | <generator> | <subword> | <commutator>
%   <subword> ::= '(' <nonEmptyWord> ')'
%   <commutator> ::= '[' <nonEmptyWord> ',' <nonEmptyWord> ']'
%   <generator> ::= corresponds to the regexp [A-Za-z][A-Za-z0-9_]*

% The commutator $[x,y]$ can be defined in two ways:
%
% - $[x, y] = x y x^{-1} y^{-1}$, in which case set `commutatorStartsWithInverses` to false
% - $[x, y] = x^{-1} y^{-1} x y$, in which case set `commutatorStartsWithInverses` to true
%
% Generator names can contain alphanumerical and punctuations characters, with the exception of
% ``^()[,]`` which are reserved.
%
% Example:
%   >>> W = replab.fp.Parser({'x' 'y'}, true);
%   >>> [ok, word] = W.parse('[x,y]^2');
%   >>> word
%       [-1, -2, 1, 2, -1, -2, 1, 2]
%
% Note:
%
%   Before parsing, the word is split into tokens. The tokens structure is a matrix with
%   two rows and as many columns as tokens. The value on the first row gives the token type,
%   while the value on the second row given the token attached data.
%
%   We have the following types:
%
%   1. <exponent>, with the data giving the value of the exponent (EXPONENT)
%   2. opening parenthesis '(' (OPEN_PAREN)
%   3. closing parenthesis ')' (CLOSE_PAREN)
%   4. opening bracket '[' (OPEN_BRACKET)
%   5. comma ',' (COMMA)
%   6. closing bracket ']' (CLOSE_BRACKET)
%   7. equality sign (EQUALITY)
%   8. identity '1' (IDENTITY)
%   9. multiplication operator '*' (TIMES)
%   10. division operator '/' (DIVIDE)
%   11. generator, with the data giving the index of the generator
%   12. end of word

    properties
        commutatorStartsWithInverses % (logical): If true, ``[x, y] = x^-1 y^-1 x y``, if false ``[x, y] = x y x^-1 y^-1``
        types % (struct): Constant values for token types
        specialChars % (charstring): List of special characters for the lexer
    end

    methods

        function self = Parser(commutatorStartsWithInverses)
            if nargin < 1
                commutatorStartsWithInverses = true;
            end
            self.commutatorStartsWithInverses = commutatorStartsWithInverses;
            self.specialChars = '^()[,]=1*/';
            self.types = struct('EXPONENT', 1, 'OPEN_PAREN', 2, 'CLOSE_PAREN', 3, ...
                                'OPEN_BRACKET', 4, 'COMMA', 5, 'CLOSE_BRACKET', 6, ...
                                'EQUALITY', 7, 'IDENTITY', 8, ...
                                'TIMES', 9, 'DIVIDE', 10, 'GENERATOR', 11, 'END', 12);
        end

        function [pos, w] = part(self, tokens, pos)
        % Parses a part
            w = [];
            % parse the part
            type = tokens(1, pos);
            types = self.types;
            if type == types.IDENTITY % identity
                pos = pos + 1;
                % w = [] still empty
            elseif type == types.GENERATOR % generator
                w = tokens(2, pos);
                pos = pos + 1;
            elseif type == types.OPEN_PAREN % subword
                pos = pos + 1;
                [pos, w] = self.nonEmptyWord(tokens, pos);
                if pos == 0 || tokens(1, pos) ~= types.CLOSE_PAREN
                    return
                end
                pos = pos + 1;
            elseif type == types.OPEN_BRACKET % commutator
                pos = pos + 1;
                [pos, lhs] = self.nonEmptyWord(tokens, pos);
                if pos == 0 || tokens(1, pos) ~= types.COMMA
                    return
                end
                pos = pos + 1;
                [pos, rhs] = self.nonEmptyWord(tokens, pos);
                if pos == 0 || tokens(1, pos) ~= types.CLOSE_BRACKET
                    return
                end
                pos = pos + 1;
                lr = [lhs rhs];
                lIrI = fliplr(-[rhs lhs]);
                if self.commutatorStartsWithInverses
                    w = [lIrI lr];
                else
                    w = [lr lIrI];
                end
            else
                pos = 0;
            end
        end

        function [pos, w] = component(self, tokens, pos)
        % Parses a component, which is a part with an optional exponent
            [pos, w] = self.part(tokens, pos);
            if pos == 0
                return
            end
            types = self.types;
            if tokens(1, pos) == types.EXPONENT
                e = tokens(2, pos);
                pos = pos + 1;
            else
                e = 1;
            end
            if e < 0
                e = -e; % invert
                w = fliplr(-w);
            end
            if e > 1
                w = repmat(w, 1, e);
            end
        end

        function [pos, w] = nonEmptyWord(self, tokens, pos)
        % Parses a word from tokens, which is a sequence of components
            [pos, w] = self.component(tokens, pos);
            if pos == 0
                return
            end
            types = self.types;
            while 1
                if tokens(1, pos) == types.TIMES
                    pos = pos + 1;
                    [pos, next] = self.component(tokens, pos);
                    if pos == 0 % times MUST be followed by component
                        w = [];
                        return
                    end
                    w = [w next];
                elseif tokens(1, pos) == types.DIVIDE
                    pos = pos + 1;
                    [pos, next] = self.component(tokens, pos);
                    if pos == 0 % divide MUST be followed by component
                        w = [];
                        return
                    end
                    w = [w -fliplr(next)];
                else
                    [nextPos, next] = self.component(tokens, pos);
                    if nextPos == 0
                        % success, we're at the end of our empty word
                        return
                    end
                    pos = nextPos;
                    w = [w next];
                end
            end
        end

        function [pos, w] = word(self, tokens, pos)
        % Parses a word
            w = [];
            [pos1, w1] = self.nonEmptyWord(tokens, pos);
            if pos1 == 0
                return % pos stays in place
            else
                pos = pos1;
                w = w1;
            end
        end

        function [pos, relators] = equation(self, tokens, pos)
        % Parses an equation and returns the relators as a cell array
        %
            r = [];
            relators = [];
            types = self.types;
            [pos, lhs] = self.nonEmptyWord(tokens, pos);
            if pos == 0
                return
            end
            if tokens(1, pos) ~= types.EQUALITY
                % this is already a relator, return
                relators = {lhs};
            end
            elements = {lhs};
            while tokens(1, pos) == types.EQUALITY
                pos = pos + 1;
                [pos, w] = self.nonEmptyWord(tokens, pos);
                if pos == 0
                    return
                end
                elements{1,end+1} = w;
            end
            rhs = -fliplr(elements{end});
            relators = cellfun(@(l) [l rhs], elements(1:end-1), 'uniform', false);
        end

        function [pos, relators] = relations(self, tokens, pos)
        % Parses a non-empty sequence of relations separated by commas
        %
        % Returns the relators as a cell array
            res = {};
            relators = [];
            types = self.types;
            [pos r] = self.equation(tokens, pos);
            if pos == 0
                return
            end
            res = r;
            while tokens(1, pos) == types.COMMA
                pos = pos + 1;
                [pos r] = self.equation(tokens, pos);
                if pos == 0
                    return
                end
                res = horzcat(res, r);
            end
            relators = res;
        end

        function [ok, names, relatorLetters] = parsePresentation(self, str)
        % Parses a string representing a group presentation
        %
        % Args:
        %   str (charstring): String to parse
        %
        % Returns
        % -------
        %   ok:
        %     logical: Whether the parse was successful
        %   names:
        %     cell(1,\*) of charstring: Names of the generators
        %   relatorLetters:
        %     cell(1,\*) of integer(1,\*): Relators given as the letters composing words
            ok = false;
            names = [];
            types = self.types;
            relatorLetters = [];
            str = strtrim(str);
            if str(1) ~= '<' || str(end) ~= '>'
                return
            end
            str = str(2:end-1); % cut the < >
            parts = strsplit(str, '|');
            nameStr = parts{1};
            names = cellfun(@strtrim, strsplit(strtrim(nameStr), {' ', ','}), 'uniform', false);
            if length(parts) < 2
                relStr = '';
            else
                relStr = parts{2};
            end
            relStr = strtrim(relStr);
            if isempty(relStr)
                relatorLetters = {};
            else
                [ok, tokens] = self.lex(relStr, names);
                if ~ok
                    names = [];
                    return
                end
                pos = 1;
                [pos, relatorLetters] = self.relations(tokens, pos);
                if pos == 0 || tokens(1, pos) ~= types.END
                    names = [];
                    relators = [];
                    return
                end
            end
            ok = true;
        end

        function [ok, T] = lex(self, str, names)
        % Splits the given string into tokens
        %
        % See the class documentation string for `.Parser`
            ok = false; % if we return midway, it's because of an invalid string, so set ok = false by default
            T = zeros(2, 0);
            n = length(str);

            % whitespace mask
            ws = [isstrprop(str, 'wspace') false];

            types = self.types;
            spec = self.specialChars;

            % valid characters in a generator name
            gstart = isstrprop(str, 'alpha');
            gnext = isstrprop(str, 'alphanum') | str == '_';
            gstart = [gstart false]; % add a false value so that the find below always returns something
            gnext = [gnext false];
            % digits
            dg = [isstrprop(str, 'digit') false];

            i = 1;
            while i <= n
                % skip whitespace
                i = find(~ws(i:end), 1) + i - 1;
                if i > n
                    break
                end
                type = find(spec == str(i), 1);
                if isempty(type)
                    if gstart(i) % is a generator
                        j = find(~gnext(i+1:end), 1) + i; % find the end of the name block
                        name = str(i:j-1);
                        [res, ind] = ismember(name, names);
                        if ~res
                            return
                        end
                        token = [types.GENERATOR; ind];
                        T = [T token];
                        i = j;
                    else % is not a generator, so unrecognized
                        return
                    end
                elseif type == 1 % power
                    i = i + 1;
                    i = find(~ws(i:end), 1) + i - 1; % skip whitespace
                    if i > n
                        return
                    end
                    if str(i) == '-'
                        sign = -1;
                        i = i + 1;
                    else
                        sign = 1;
                    end
                    if i > n || ~dg(i)
                        return
                    end
                    j = find(~dg(i+1:end), 1) + i; % find the end of the integer
                    e = str2num(str(i:j-1)) * sign;
                    token = [type; e];
                    T = [T token];
                    i = j;
                else
                    token = [type; 0];
                    T = [T token];
                    i = i + 1;
                end
            end
            token = [types.END; 0];
            T = [T token];
            ok = true;
        end

    end

end
