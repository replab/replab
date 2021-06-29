classdef Parser < replab.Str
% Parses a cyclotomic expression
%
% For now, it returns a floating-point approximation of the represented cyclotomic number.
%
% Here is the grammar of the expressions it parses:
%
% ::
%
%    <expression> ::= <sum> <end>
%    <sum> ::= <product> ( ('+' <product>) | ('-' <product>) )*
%    <product> ::= <factor> ( ('*' <factor>) | ('/' <factor>) )*
%    <factor> ::= <atom> | ('-' <factor>) | ('+' <factor>) | (<atom> '^' '-'? <factor>)
%    <atom> ::= (<integer>) | ( 'E' '(' '-'? <integer> ')' ) | ( '(' <sum> ')' )
%    <integer> ::= '0' | ( ['1'-'9'] ['0'-'9']* )
%
% Example:
%   >>> replab.cyclo.Parser.parse('2+2')
%       4
%
% Note:
%
%   Before parsing, the expression is split into tokens. The tokens structure is a matrix with two rows and as many columns
%   as tokens. The value on the first row gives the token type, while the value on the second row given the token attached data,
%   if necessary.

    properties (Constant = true)
        types = struct('PLUS', 1, 'MINUS', 2, 'TIMES', 3, 'DIVIDE', 4, 'POWER', 5, 'OPEN_PAREN', 6, 'CLOSE_PAREN', 7, 'E', 8, 'INTEGER', 9, 'END', 10);
        specialChars = '+-*/^()E';
    end

    methods (Static)

        function value = parse(str)
            value = [];
            [ok, tokens] = replab.cyclo.Parser.lex(str);
            if ~ok
                error('Invalid expression %s', str);
            end
            [pos, value] = replab.cyclo.Parser.parseExpression(tokens, 1);
            if pos == 0
                error('Invalid expression %s', str);
            end
        end

        function [pos, value] = parseExpression(tokens, pos)
            value = [];
            [pos, value] = replab.cyclo.Parser.parseSum(tokens, pos);
            if ~pos
                return
            end
            if tokens(1, pos) ~= replab.cyclo.Parser.types.END
                pos = 0;
                return
            end
        end


        function [pos, value] = parseSum(tokens, pos)
            value = [];
            types = replab.cyclo.Parser.types;
            [pos, lhs] = replab.cyclo.Parser.parseProduct(tokens, pos);
            if pos == 0
                return
            end
            while tokens(1, pos) == types.PLUS || tokens(1, pos) == types.MINUS
                if tokens(1, pos) == types.PLUS
                    s = 1;
                else
                    s = -1;
                end
                pos = pos + 1;
                [pos, rhs] = replab.cyclo.Parser.parseProduct(tokens, pos);
                if pos == 0
                    return
                end
                lhs = lhs + s * rhs;
            end
            value = lhs;
        end

        function [pos, value] = parseProduct(tokens, pos)
            value = [];
            types = replab.cyclo.Parser.types;
            [pos, lhs] = replab.cyclo.Parser.parseFactor(tokens, pos);
            if pos == 0
                return
            end
            while tokens(1, pos) == types.TIMES || tokens(1, pos) == types.DIVIDE
                isDivide = (tokens(1, pos) == types.DIVIDE);
                pos =  pos + 1;
                [pos, rhs] = replab.cyclo.Parser.parseFactor(tokens, pos);
                if pos == 0
                    return
                end
                if isDivide
                    lhs = lhs / rhs;
                else
                    lhs = lhs * rhs;
                end
            end
            value = lhs;
        end

        function [pos, value] = parseFactor(tokens, pos)
            types = replab.cyclo.Parser.types;
            value = [];
            if tokens(1, pos) == types.PLUS
                pos = pos + 1;
                [pos value] = replab.cyclo.Parser.parseFactor(tokens, pos);
                if pos == 0
                    value = [];
                end
            elseif tokens(1, pos) == types.MINUS
                pos = pos + 1;
                [pos value] = replab.cyclo.Parser.parseFactor(tokens, pos);
                value = -value;
                if pos == 0
                    value = [];
                end
            else
                [pos value] = replab.cyclo.Parser.parseAtom(tokens, pos);
                if pos == 0
                    value = []
                    return
                end
                if tokens(1, pos) == types.POWER
                    pos = pos + 1;
                    if tokens(1, pos) == types.MINUS
                        pos = pos + 1;
                        s = -1;
                    else
                        s = 1;
                    end
                    if tokens(1, pos) == types.INTEGER
                        value = value^(s*tokens(2, pos));
                        pos = pos + 1;
                    else
                        pos = [];
                        value = [];
                        return
                    end
                end
            end
        end

        function [pos value] = parseE(tokens, pos)
            value = [];
            types = replab.cyclo.Parser.types;
            if tokens(1, pos) ~= types.E
                pos = 0;
                return
            end
            pos = pos + 1;
            if tokens(1, pos) ~= types.OPEN_PAREN
                pos = 0;
                return
            end
            pos = pos + 1;
            if tokens(1, pos) == types.MINUS
                s = -1;
                pos = pos + 1;
            else
                s = 1;
            end
            if tokens(1, pos) ~= types.INTEGER
                pos = 0;
                return
            end
            value = exp(2i*pi/tokens(2, pos));
            pos = pos + 1;
            if tokens(1, pos) ~= types.CLOSE_PAREN
                pos = 0;
                value = [];
                return
            end
            pos = pos + 1;
        end

        function [pos value] = parseAtom(tokens, pos)
            types = replab.cyclo.Parser.types;
            value = [];
            if tokens(1, pos) == types.INTEGER
                value = tokens(2, pos);
                pos = pos + 1;
            elseif tokens(1, pos) == types.E
                [pos value] = replab.cyclo.Parser.parseE(tokens, pos);
                if pos == 0
                    return
                end
            elseif tokens(1, pos) == types.OPEN_PAREN
                pos = pos + 1;
                [pos value] = replab.cyclo.Parser.parseSum(tokens, pos);
                if tokens(1, pos) == ~types.CLOSE_PAREN
                    pos = [];
                    value = [];
                end
            else
                pos = 0;
            end
        end

        function [ok, T] = lex(str)
        % Splits the given string into tokens
        %
        % See the class documentation string for `.Parser`
            ok = false; % if we return midway, it's because invalid string, so ok is false by default
            T = zeros(2, 0);
            n = length(str);

            % strip whitespace
            str = str(~isstrprop(str, 'wspace'));

            types = replab.cyclo.Parser.types;
            specialChars = replab.cyclo.Parser.specialChars;

            % valid integers
            dg = [isstrprop(str, 'digit') false];

            i = 1;
            while i <= n
                type = find(specialChars == str(i), 1);
                if isempty(type)
                    if dg(i) % is an integer
                        j = find(~dg(i+1:end), 1) + i; % find the end of the integer
                        intstr = str(i:j-1);
                        if intstr(1) == '0'
                            if length(intstr) > 1
                                return % '0' is the only integer that can start with the digit '0'
                            end
                            token = [types.INTEGER; 0];
                        else
                            token = [types.INTEGER; str2num(intstr)];
                        end
                    else % is not an integer, so unrecognized
                        return
                    end
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
