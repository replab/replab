classdef CyclotomicTable < replab.Str

    properties
        matrix % (`+replab.cyclotomic` (n,n)): Coefficients
        variables % (cell(1,nV) of charstring): Variables
        values % (`+replab.cyclotomic` (1,nV)): Values
        hashes % (integer(1,\*)): Hashes of the values
    end

    methods

        function self = CyclotomicTable(matrix)
            self.matrix = matrix;
            self.variables = cell(1, 0);
            self.values = replab.cyclotomic.zeros(1, 0);
            self.hashes = zeros(1, 0);
        end

        function [ind negate conjugate] = find(self, c)
            tests = [c -c conj(c) -conj(c)];
            negates = [false true false true];
            conjugates = [false false true true];
            for i = 1:length(tests)
                c = tests(i);
                h = c.hash;

            end
        end

    end

end
