classdef RealRep < replab.Str
% A finitely generated group real representation on GL_d(R)
    
    properties (SetAccess = protected)
        group; % Group represented
        dimension; % Representation dimension
        images; % Generator images
        invImages; % Generator inverse images
    end
    
    properties (Access = protected)
        M;
        centralizerAlgebra_;
    end
    
    methods
        
        function self = RealRep(group, dimension, images, invImages)
        % Constructs a representation from a group's generator images
            assert(isa(group, 'replab.FinitelyGeneratedGroup'));
            assert(length(images) == group.nGenerators);
            self.group = group;
            self.dimension = dimension;
            self.images = images;
            self.invImages = invImages;
            d = dimension;
            self.M = replab.GroupFun('Mat', @isequal, @() rand(d, d), @(x, y) x*y, eye(d), @(x) inv(x));
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
        
        function rho = image(self, g)
            word = self.group.factorization(g);
            rho = self.M.identity;
            for i = 1:length(word.indices)
                ind = word.indices(i);
                e = word.exponents(i);
                if e > 0
                    ge = self.M.composeN(self.images{ind}, e);
                else
                    ge = self.M.composeN(self.invImages{ind}, -e);
                end
                rho = self.M.compose(rho, ge);
            end
        end

        function b = isUnitary(self)
            d = self.dimension;
            tol = replab.Settings.doubleEigTol;
            for i = 1:length(self.images)
                I = self.images{i};
                if ~self.M.isIdentity(I*I')
                    b = false;
                    return
                end
            end
            b = true;                    
        end
        
        function c = centralizerAlgebra(self)
            if isempty(self.centralizerAlgebra_)
                self.centralizerAlgebra_ = replab.RealCentralizerAlgebra(self);
            end
            c = self.centralizerAlgebra_;
        end

    end
    
end
