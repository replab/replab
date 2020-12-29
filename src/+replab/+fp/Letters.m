classdef Letters

    methods (Static)

        function letters = parse(word, names)
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

        function word = print(letters, names, times)
        % Prints a word from its letters
        %
        % Args:
        %   letters (integer(1,\*)): Letters of the word as generator indices
        %   names (cell(1,\*) of charstring): Generator names
        %   times (charstring, optional): Composition operator, default value ``' '``
        %
        % Returns:
        %   charstring: Word given as an explicit string
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
            letters = replab.fp.Letters.reduce(letters);
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
                        word = sprintf('%s^%d', word, e);
                    end
                end
            end
        end

        function x = reduce(x)
        % Finds the reduced form of the word described the given letters
        %
        % It eliminates all products of the form ``[i -i]`` or ``[-i i]`` for the index ``i`` of a generator.
        %
        % Args:
        %   x (integer(1,\*)): Letters
        %
        % Returns:
        %   integer(1,\*): Letters of the reduced word
            if isempty(x)
                x = zeros(1, 0);
                return
            end
            i = find(x(1:end-1) == -x(2:end), 1);
            if ~isempty(i)
                while i < length(x)
                    if x(i) == -x(i+1)
                        x = [x(1:i-1) x(i+2:end)];
                        if i > 1
                            i = i - 1;
                        end
                    else
                        i = i + 1;
                    end
                end
            end
        end

        function z = inverse(x)
        % Inverts the word described the given letters
        %
        % Args:
        %   x (integer(1,\*)): Letters
        %
        % Returns:
        %   integer(1,\*): Letters of the word inverse
            z = -fliplr(x);
        end

        function z = compose(x, y)
        % Composes the two words described the given letters
        %
        % If both arguments are reduced, it guarantees that the result is reduced too.
        %
        % Args:
        %   x (integer(1,\*)): Letters of the left hand side word
        %   y (integer(1,\*)): Letters of the right hand side word
        %
        % Returns:
        %   integer(1,\*): Letters of the product
            xe = length(x);
            ys = 1;
            yn = length(y);
            while xe > 0 && ys <= yn && x(xe) == -y(ys)
                xe = xe - 1;
                ys = ys + 1;
            end
            z = [x(1:xe) y(ys:end)]; % if empty, will give zeros(1, 0)
        end

        function r = cyclicallyReduce(r)
            r = replab.fp.Letters.reduce(r);
            while length(r) >= 2 && r(1) == -r(end)
                r = r(2:end-1);
            end
        end

        function c = cyclicConjugates(r)
            r = replab.fp.Letters.cyclicallyReduce(r);
            n = length(r);
            C = zeros(n, n);
            for i = 1:n
                C(i,:) = [r(i:end) r(1:i-1)];
            end
            C = unique(C, 'rows');
            c = cell(1, size(C, 1));
            for i = 1:size(C, 1)
                c{i} = C(i,:);
            end
        end

        function list = unique(list)
        % Returns a list of unique words
            m = max(cellfun(@length, list)) + 1; % so we have a guard
            L = zeros(length(list), m);
            for i = 1:length(list)
                l = list{i};
                L(i,1:length(l)) = l;
            end
            L = unique(L, 'rows');
            list = cell(1, size(L, 1));
            for i = 1:size(L, 1)
                l = L(i, :);
                j = find(l == 0, 1);
                list{i} = l(1:j-1);
            end
        end

    end

end
