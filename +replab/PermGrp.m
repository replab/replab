classdef PermGrp < replab.FiniteGrp
% Represents a permutation group

    methods
        
        function d = domainSize(self) 
        % Returns d when this group acts on {1, ..., d}
            error('Not implemented');
        end
        
        function b = contains(self, permutation) 
        % Returns whether "permutation" is member of this group
            error('Not implemented');
        end

        function o = order(self)
        % Returns the order of the group as a vpi value
            error('Not implemented');
        end

        function G = generatorsAsMatrix(self)     
        % Returns a double matrix whose rows are generators
            error('Not implemented');
        end
        
        function p = uniformlyRandomElement(self) 
        % Returns a random element from the group, guaranteed
        % to be distributed uniformly
            error('Not implemented');
        end
        
        function E = element(self, i) 
        % Returns the i-th element in the lexicographic order
        % i is either a double or a vpi
            error('Not implemented');
        end

        function E = elements(self) 
        % Enumerates the group permutations, sorted, in a row cell vector
            error('Not implemented');
        end
        
        function E = elementsAsMatrix(self)    
        % Returns all permutations from the group as row
        % vectors in a double matrix
            error('Not implemented');
        end
        
        function w = factorization(self, permutation) 
        % Return "permutation" as a word in the generators of this group
            error('Not implemented');
        end
        
        function p = evaluateWord(self, word); 
        % Evalutes the given word using this group generators
            error('Not implemented');
        end
        
    end
    
    methods (Static) % Factory methods
        
        function G = trivial(domainSize)
        % Constructs the trivial group on the domain
            G = replab.prv.PermGrpList.fromGeneratorsAsMatrix(zeros(0, domainSize));
        end
        
        function G = symmetric(n)
        % Constructs the symmetric group acting on n elements
            if n == 1
                G = replab.PermGrp.trivial(n);
            elseif n == 2
                G = replab.PermGrp.fromGenerators([2 1]);
            else
                G = replab.PermGrp.fromGenerators([2:n 1], [2 1 3:n]);
            end
        end
        
        function G = cyclic(n)
        % Constructs the cyclic group acting on n elements
            if n == 1
                G = replab.PermGrp.trivial(n);
            else
                G = replab.PermGrp.fromGenerators([2:n 1]);
            end
        end
        
        function G = fromGeneratorsAsMatrix(genMat)
            errmsg = 'To create the trivial group, pass zeros(0, n) so that the domain size n is known';
            assert(size(genMat, 2) > 0, errmsg);
            G = replab.prv.PermGrpList.fromGeneratorsAsMatrix(genMat);
        end
        
        function G = fromDomainSizeAndGenerators(domainSize, varargin)
            nG = length(varargin);
            genMat = zeros(nG, domainSize);
            for i = 1:nG
                genMat(i, :) = varargin{i};
            end
            G = replab.prv.PermGrpList.fromGeneratorsAsMatrix(genMat);
        end
        
        function G = fromGenerators(varargin)
            errmsg = 'Argument list cannot be empty. Use fromDomainSizeAndGenerators instead';
            assert(nargin > 0, errmsg);
            g = varargin{1};
            domainSize = length(g);
            G = replab.PermGrp.fromDomainSizeAndGenerators(domainSize, varargin{:});
        end
        
    end
    
end
