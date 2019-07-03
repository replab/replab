classdef WreathProductGroup < replab.SemidirectProductGroup
    
    properties (SetAccess = protected)
        A; % factor of base group
        n; % number of copies of the base group
    end
        
    methods
        
        function self = WreathProductGroup(H, A)
            assert(isa(H, 'replab.PermutationGroup'));
            assert(isa(A, 'replab.FiniteGroup'));
            n = H.domainSize;
            base = replab.DirectProductGroup.power(A, n);
            phi = replab.perm.PermutationCellAction(H, base);
            self = self@replab.SemidirectProductGroup(phi);
            self.n = n;
            self.A = A;
        end
        
% $$$         function phi = imprimitiveMorphism(self, phiA)
% $$$             if nargin < 2
% $$$                 phi = @(w) self.imprimitivePermutation(w);
% $$$             else
% $$$                 phi = @(w) self.imprimitivePermutation(w, phiA);
% $$$             end
% $$$         end
% $$$ 
% $$$         function p = imprimitivePermutation(self, w, phiA)
% $$$         %      w: wreath element to compute the image of
% $$$         %    phi: morphism from elements of A to permutations
% $$$         %         optional, if omitted default to identity
% $$$         %         which is valid only when A is a permutation group
% $$$             if nargin < 3
% $$$                 assert(isa(self.A, 'replab.PermutationGroup'));
% $$$                 phiA = @(x) x;
% $$$             end
% $$$             n = self.n;
% $$$             h = w{1};
% $$$             base = w{2};
% $$$             im = phiA(base{1});
% $$$             d = length(im);
% $$$             basePerm = im;
% $$$             shift = d;
% $$$             for i = 2:n
% $$$                 im = phiA(base{i}) + shift;
% $$$                 basePerm = [basePerm im];
% $$$                 shift = shift + d;
% $$$             end
% $$$             p = reshape(basePerm, [d n]);
% $$$             p(:,h) = p;
% $$$             p = p(:)';
% $$$         end
% $$$         
% $$$         function p = primitiveMorphism(self, phi)
% $$$         end
% $$$         
% $$$         function p = primitiveImage(self, w, phi)
% $$$         end
        
        function rep = imprimitiveRep(self, Arep)
            rep = replab.rep1.WreathProductImprimitiveRep(self, Arep);
        end
        
        function rep = primitiveRep(self, Arep)
            rep = replab.rep1.WreathProductPrimitiveRep(self, Arep);
        end
        
    end
    
end
