classdef IrreducibleRep < replab.Rep
% Expresses the block-diagonal representation corresponding to an irreducible decomposition
%
% It uses the known structure to provide a cleaner expression of the represenation.
   
    properties
        irreducible % replab.Irreducible: Decomposition of the representation into irreducibles
    end
    
    methods
        
        function self = IrreducibleRep(irreducible)
        % Constructs a representation from an irreducible decomposition
            parent = irreducible.parent;
            assert(isequal(parent.isUnitary, true));
            self = self@replab.Rep(parent.group, parent.field, parent.dimension);
            self.irreducible = irreducible;
        end
        
        function rho = image(self, g)
            I = self.irreducible;
            % Construct the blocks in the block diagonal image
            blocks = {};
            for i = 1:I.nComponents
                % For each isotypic component
                C = I.component(i);
                % Average over all copies to reduce noise
                Cj = C.copy(1);
                block = Cj.image(g);
                for j = 2:C.multiplicity
                    Cj = C.copy(j);
                    block = block + Cj.image(g);
                end
                block = block/C.multiplicity;
                % And add the correct number of duplicates 
                for j = 1:C.multiplicity
                    blocks{1,end+1} = block;
                end
            end
            % Computes the block diagonal basis
            rho = blkdiag(blocks{:});
        end
        
    end
    
end
