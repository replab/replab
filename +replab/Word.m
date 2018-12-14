classdef Word
% An associative word of the form 
% w = x(i1)^e1 x(i2)^e2 ...
% where indices = [i1 i2 ...] and exponents = [e1 e2 ...]
%
% Exponents are nonzero signed integers
    
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
        
        function l = length(self)
        % Returns the length of the word, where each
        % product by a generator or its inverse counts as one letter
            l = sum(abs(self.exponents));
        end
        
        function b = isGenerator(self)
            b = (length(self.indices) == 1) && (self.exponents == 1);
        end
        
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
            W = replab.Word.fromIndicesAndExponents(newInd, newExp);
        end
        
        function W = mrdivide(self, rhs)
            W = self * rhs.inverse;
        end
        
        function W = inv(self)
            newInd = fliplr(self.indices);
            newExp = fliplr(-self.exponents);
            W = replab.Word(newInd, newExp);
        end
        
        function y = mpower(self, e)
            if e < 0
                y = self.mpower(self.inv, -e);
            elseif e == 0
                y = replab.Word.identity;
            else
                x = self;
                y = replab.Word.identity;
                while e > 1
                    if mod(e, 2) == 0 % n even
                        x = x * x;
                        e = e / 2;
                    else
                        y = x * y;
                        x = x * x;
                        e = (e - 1)/2;
                    end
                end
                y = x * y;
            end
        end
        
    end
        
    methods (Static) % CONSTRUCTION METHODS
        
        function W = identity
            W = replab.Word([], []);
        end
        
        function W = generator(index)
            W = replab.Word(index, 1);
        end
        
        function W = random(nGenerators, maxLength)
            l = randi(maxLength + 1) - 1;
            % exponents are +1 or -1
            e = randi(2, 1, l);
            e(e == 2) = -1;
            g = randi(nGenerators, 1, l); % generator indices
            W = replab.Word.fromIndicesAndExponents(g, e);
        end
        
        function W = fromIndicesAndExponents(indices, exponents)
        % Reduces the given indices and exponents
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
