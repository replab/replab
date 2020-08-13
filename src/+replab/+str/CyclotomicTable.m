classdef CyclotomicTable < replab.Str

    properties
        matrix % (`+replab.cyclotomic` (m,n)): Coefficients
        strings % (cell(m,n) of charstring): String representation of coefficients
        variables % (cell(1,nV) of charstring): Variables
        values % (`+replab.cyclotomic` (1,nV)): Values
        hashes % (integer(1,\*)): Hashes of the values
    end

    methods (Static)

        function s = variableName(n)
            if n <= 26
                s = char('A' + n - 1);
            else
                s = sprintf('_%d', n);
            end
        end

    end


    methods

        function self = CyclotomicTable(matrix)
            self.matrix = matrix;
            self.variables = cell(1, 0);
            self.values = replab.cyclotomic.zeros(1, 0);
            self.hashes = zeros(1, 0);
            self.strings = cell(size(matrix));
            for i = 1:size(matrix, 1)
                for j = 1:size(matrix, 2)
                    c = matrix(i,j);
                    if c == 1i
                        self.strings{i,j} = 'i';
                    elseif c == -1i
                        self.strings{i,j} = '-i';
                    elseif c.isWhole
                        self.strings{i,j} = sprintf('%d', double(c));
                    else
                        [ind negate conjugate] = self.find(matrix(i, j));
                        s = self.variables{ind};
                        if negate
                            s = ['-' s];
                        end
                        if conjugate
                            s = [s '*'];
                        else
                            s = [s ' '];
                        end
                        self.strings{i,j} = s;
                    end
                end
            end
        end

        function [ind negate conjugate] = find(self, c)
            for negate = [false true]
                c1 = c;
                if negate
                    c1 = -c1;
                end
                for conjugate = [false true]
                    c2 = c1;
                    if conjugate
                        c2 = conj(c2);
                    end
                    h = c2.hash;
                    ind = find(self.hashes == h);
                    ind = ind(find(self.values(ind) == c2));
                    if ~isempty(ind)
                        return
                    end
                end
            end
            ind = length(self.variables) + 1;
            self.variables = [self.variables replab.str.CyclotomicTable.variableName(ind)];
            self.values = [self.values c];
            self.hashes = [self.hashes c.hash];
            negate = false;
            conjugate = false;
        end

    end

end
