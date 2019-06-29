classdef RealIrrep < replab.RealRep
    
    properties (SetAccess = protected)
        images1;    % Images of a single copy 
        imagesInv1; % 
        dimension1; % Dimension of a single copy of this irreducible representation
        multiplicity; % Representation multiplicity
        divisionAlgebra; % Division algebra (R, C or H) as the type of this real irreducible rep.
    end
    
    methods
        
        function self = RealIrrep(group, dimension1, multiplicity, divisionAlgebra, images1, imagesInv1, parent, U, Uinv)
            if nargin < 7
                parent = [];
                U = [];
                Uinv = [];
            end
            dimension = dimension1 * multiplicity;
            images = cellfun(@(X) kron(eye(multiplicity), X), images1, 'UniformOutput', false);
            imagesInv = cellfun(@(X) kron(eye(multiplicity), X), imagesInv1, 'UniformOutput', false);
            self = self@replab.RealRep(group, dimension, images, imagesInv, parent, U, Uinv);
            self.images1 = images1;
            self.imagesInv1 = imagesInv1;
            self.dimension1 = dimension1;
            self.multiplicity = multiplicity;
            self.divisionAlgebra = divisionAlgebra;
        end
        
        function sub = copy(self, m)
        % Returns the m-th copy of this representation as an irreducible representation of
        % multiplicity 1
            if isequal(self.parent, [])
                parent1 = [];
                U1 = [];
                Uinv1 = [];
            else
                parent1 = self.parent;
                ind = (m-1)*self.dimension1 + (1:self.dimension1);
                U1 = self.U(:, ind);
                Uinv1 = self.Uinv(ind, :);
            end
            sub = replab.RealIrrep(self.group, self.dimension1, 1, self.divisionAlgebra, ...
                                   self.images1, self.imagesInv1, ...
                                   parent1, U1, Uinv1);
        end
        
        function s = headerStr(self)
            switch self.divisionAlgebra.shortName
              case 'R'
                adj = 'real-type';
              case 'C'
                adj = 'complex-type';
              case 'H'
                adj = 'quaternionic-type';
            end
            if self.isUnitary
                t = sprintf('Unitary %s irreducible representation', adj);
            else
                t = sprintf('Irreducible %s representation', adj);
            end
            s = sprintf('%s of multiplicity %d, dimension %d', t, self.multiplicity, self.dimension1);
        end
        
        function c = centralizerAlgebra(self)
            c = replab.RealIrreducibleCentralizerAlgebra(self);
        end
        
    end

    methods (Static)
        
        function rep = fromParentRealRep(parent, dimension1, multiplicity, divisionAlgebra, U, Uinv)
        % Constructs a subrepresentation of (single copy) dimension "dimension1" and given "multiplicity" 
        % in the given division algebra from a "parent" representation and change of basis matrices U and Uinv, 
        % such that subrep.image(g) = Uinv*parent.image(g)*U
            nG = parent.group.nGenerators;
            images1 = cell(1, nG);
            imagesInv1 = cell(1, nG);
            for i = 1:nG
                rho = parent.images{i};
                rhoInv = parent.imagesInv{i};
                img = zeros(dimension1, dimension1);
                imgInv = zeros(dimension1, dimension1);
                for j = 1:multiplicity
                    ind = (j-1)*dimension1 + (1:dimension1);
                    img = img + Uinv(ind,:)*rho*U(:,ind);
                    imgInv = imgInv + Uinv(ind,:)*rhoInv*U(:,ind);
                end
                img = img / multiplicity;
                imgInv = imgInv / multiplicity;
                img = divisionAlgebra.projectMatrix(img);
                imgInv = divisionAlgebra.projectMatrix(imgInv);
                images1{i} = img;
                imagesInv1{i} = imgInv;
            end
            rep = replab.RealIrrep(parent.group, dimension1, multiplicity, divisionAlgebra, ...
                                   images1, imagesInv1, parent, U, Uinv);
        end
        
    end
    
end
