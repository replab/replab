classdef RealIrrep < replab.RealRep
    properties (SetAccess = protected)
        images1;
        imagesInv1;
        dimension1;
        multiplicity;
        divisionAlgebraSize;
        divisionAlgebra;        
    end
    methods
        function self = RealIrrep(group, dimension1, multiplicity, divisionAlgebra, images1, imagesInv1)
            self.images = cellfun(@(X) kron(eye(multiplicity), X), images1, 'UniformOutput', false);
            self.imagesInv = cellfun(@(X) kron(eye(multiplicity), X), imagesInv1, 'UniformOutput', false);
            self.group = group;
            self.dimension1 = 
        end
        function sub = single(self)
            
        end
        function s = str(self)
            if self.isUnitary
                t = 'Unitary representation';
            else
                t = 'Representation';
            end
            s = sprintf('%s of dimension %d with generator images', t, self.dimension);
            for i = 1:length(self.images)
                gen = char('a' + i - 1);
                s = [s char(10) '- ' gen ':' char(10)];
                img = replab.prependLines(replab.strOf(self.images{i}), '    ');
                s = [s img char(10)];
            end
        end
        function b = isUnitary(self)
        % Specialize
            b = all(cellfun(@(U) self.M.isIdentity(U*U'), self.images1));
        end
        
        function b = isMonomial(self)
            b = all(cellfun(@(U) replab.SignedPermutations.isSignedPermutationMatrix(U), self.images1));
        end
        function 
    end
end
