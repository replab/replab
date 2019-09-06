classdef Word < replab.Str
% An associative word of the form 
% w = x(i1)^e1 x(i2)^e2 ...
% where indices = [i1 i2 ...] and exponents = [e1 e2 ...]
%
% Exponents are nonzero signed integers
    
    properties (SetAccess = protected)
        indices; % 1 x L double
        exponents; % 1 x L double
    end
    
    methods (Access = private)
        
        function self = Word(indices, exponents)
            self.indices = indices;
            self.exponents = exponents;
        end
        
    end
    
    methods (Static)

        function b = eqv(w1, w2)
            b = isequal(w1.indices, w2.indices) && isequal(w1.exponents, w2.exponents);
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
        
        function s = headerStr(self)
            s = ['Word ''' self.shortStr(replab.Settings.maxColumns) ''''];
        end
        
        function s = shortStr(self, maxColumns)
            sep = '';
            s = [];
            for i = 1:length(self.indices)
                letter = char('a' + self.indices(i) - 1);
                e = self.exponents(i);
                if e == 1
                    expString = '';
                else
                    expString = ['^' num2str(e)];
                end
                s = [s sep letter expString];
                sep = ' ';
            end
        end
        
    end
    
    methods
        
        function W = mtimes(self, rhs)
            lhsInd = self.indices;
            lhsExp = self.exponents;
            rhsInd = rhs.indices;
            rhsExp = rhs.exponents;
            nL = length(lhsInd);
            nR = length(rhsInd);
            l = length(lhsInd);
            r = 1;
            while true
                if l == 0 || r > nR || lhsInd(l) ~= rhsInd(r)
                    midInd = [];
                    midExp = [];
                    break
                end
                if lhsExp(l) ~= -rhsExp(r)
                    midInd = lhsInd(l);
                    midExp = lhsExp(l) + rhsExp(r);
                    l = l - 1;
                    r = r + 1;
                    break
                end
                l = l - 1;
                r = r + 1;
            end
            newInd = [lhsInd(1:l) midInd rhsInd(r:end)];
            newExp = [lhsExp(1:l) midExp rhsExp(r:end)];
            W = replab.Word(newInd, newExp);
        end
        
        function W = mrdivide(self, rhs)
            W = self * inv(rhs);
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
        
        function W = generator(i)
            W = replab.Word(i, 1);
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
                    if i > 1
                        i = i - 1;
                    end
                elseif exponents(i) == 0
                    indices = [indices(:,1:i-1) indices(:,i+1:end)];
                    exponents = [exponents(:,1:i-1) exponents(:,i+1:end)];
                else
                    i = i + 1;
                end
            end
            W = replab.Word(indices, exponents);
        end
        
    end
    
end
