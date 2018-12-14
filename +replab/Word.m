classdef Word
    
    properties
        indices;
        exponents;
    end
    
    methods (Access = private)
        
        function self = Word(indices, exponents)
            self.indices = indices;
            self.exponents = exponents;
        end
        
    end
    
    methods
        
        function disp(self)
            disp(self.toString);
        end
        
        function str = toString(self)
            sep = '';
            str = [];
            for i = 1:length(self.indices)
                letter = char('a' + self.indices(i) - 1);
                e = self.exponents(i);
                if e == 1
                    expString = '';
                else
                    expString = ['^' num2str(e)];
                end
                str = [str sep letter expString];
                sep = ' ';
            end
        end
        
    end
    
    methods
        
        function W = mtimes(self, rhs)
            if isnumeric(rhs)
                if rhs > 0
                    rhsInd = rhs;
                    rhsExp = 1;
                else
                    rhsInd = -rhs;
                    rhsExp = -1;
                end
            else
                rhsInd = rhs.indices;
                rhsExp = rhs.exponents;
            end
            newInd = [self.indices rhsInd];
            newExp = [self.exponents rhsExp];
            W = replab.Word.reduce(newInd, newExp);
        end
        
    end
        
    methods (Static) % CONSTRUCTION METHODS
        
        function W = identity
            W = replab.Word([], []);
        end
        
        function W = generator(index, exponent)
            if nargin < 2
                exponent = 1;
            end
            if exponent == 0
                W = replab.Word.identity;
            else
                W = replab.Word(index, exponent);
            end
        end
        
        function W = random(nGenerators, maxLength)
            l = randi(maxLength + 1) - 1;
            % exponents are +1 or -1
            e = randi(2, 1, l);
            e(e == 2) = -1;
            g = randi(nGenerators, 1, l); % generator indices
            W = replab.Word.reduce(g, e);
        end
        
        function W = reduce(indices, exponents)
            i = 1;
            while i < length(indices);
                if indices(i) == indices(i+1)
                    newE = exponents(i) + exponents(i+1);
                    if newE == 0
                        newI = [];
                        newE = [];
                    else
                        newI = indices(i);
                    end
                    indices = [indices(:,1:i-1) newI indices(:,i+2:end)];
                    exponents = [exponents(:,1:i-1) newE exponents(:,i+2:end)];
                else
                    i = i + 1;
                end
            end
            W = replab.Word(indices, exponents);
        end
        
    end
    
end
