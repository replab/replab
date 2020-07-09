classdef Word < replab.Str

    properties
        group % (`+replab.FPGroup`): Group this word is part of
        letters % (integer(1,\*)): Word contents
    end

    methods

        function self = Word(group, letters)
            self.group = group;
            self.letters = replab.Word.reduceLetters(letters);
        end

        function l = length(self)
            l = length(self.letters);
        end

        function s = word2str(self)
            if self.length == 0
                s = '1';
                return
            end
            s = '';
            sep = '';
            i = 1;
            word = self.letters;
            names = self.group.names;
            while i <= length(word)
                e = sign(word(i)); % exponent
                l = abs(word(i));
                assert(e ~= 0);
                i = i + 1;
                while i <= length(word) && abs(word(i)) == l
                    e = e + sign(word(i));
                    i = i + 1;
                end
                if e ~= 0
                    s = [s sep names{l}];
                    sep = ' ';
                    if e ~= 1
                        s = [s '^' num2str(e)];
                    end
                end
            end
        end

        function s = headerStr(self)
            s = self.word2str;
        end

        function s = shortStr(self, maxColumns)
            s = self.word2str;
        end

        function s = longStr(self, maxRows, maxColumns)
            s = {self.word2str};
        end

        function l = eq(self, rhs)
            l = false;
            if self.group.groupId ~= rhs.group.groupId
                return
            end
            l = self.group.eqv(self, rhs);
        end

        function z = mtimes(self, rhs)
            assert(self.group.groupId == rhs.group.groupId, 'Must multiply words of the same group');
            z = replab.Word(self.group, self.composeLetters(self.letters, rhs.letters));
        end

        function z = inv(self)
            z = replab.Word(self.group, self.inverseLetters(self.letters));
        end

        function z = mrdivide(self, rhs)
            z = self * inv(rhs);
        end

        function z = mpower(self, n)
            z = self.group.composeN(self, n);
        end

        function g = computeImage(self, target, generatorImages)
            letters = self.letters;
            g = target.identity;
            for i = 1:length(letters)
                l = letters(i);
                if l > 0
                    g = target.compose(g, generatorImages{l});
                else
                    g = target.composeWithInverse(g, generatorImages{-l});
                end
            end
        end

    end

    methods (Static)

        function z = inverseLetters(x)
            z = -fliplr(x);
        end

        function z = composeLetters(x, y)
            xe = length(x);
            ys = 1;
            yn = length(y);
            while xe > 0 && ys <= yn && x(xe) == -y(ys)
                xe = xe - 1;
                ys = ys + 1;
            end
            z = [x(1:xe) y(ys:end)];
        end

        function x = reduceLetters(x)
            i = find(x(1:end-1) == -x(2:end));
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

    end

end
