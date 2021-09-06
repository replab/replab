classdef SemidirectProductGroup_finite < replab.SemidirectProductGroup & replab.NiceFiniteGroup

    methods

        function self = SemidirectProductGroup_finite(phi, type, varargin)
            assert(isa(phi, 'replab.Action'));
            H = phi.G;
            N = phi.P;
            assert(isa(N, 'replab.FiniteGroup'));
            assert(isa(H, 'replab.FiniteGroup'));
            identity = {H.identity N.identity};
            generators = cell(1, H.nGenerators + N.nGenerators);
            for i = 1:length(H.generators)
                generators{i} = {H.generator(i) N.identity};
            end
            for i = 1:length(N.generators)
                generators{H.nGenerators + i} = {H.identity N.generator(i)};
            end
            self@replab.NiceFiniteGroup(identity, generators, type, varargin{:});
            self.phi = phi;
            self.H = H;
            self.N = N;
        end

    end

    methods (Access = protected)

        function g = atFun(self, ind)
            ind = ind - 1;
            indN = mod(ind, self.N.order);
            indH = (ind - indN)/self.N.order;
            g = {self.H.elementsSequence.at(indH + 1) self.N.elementsSequence.at(indN + 1)};
        end

        function ind = findFun(self, g)
            indH = self.H.elementsSequence.find(g{1});
            indN = self.N.elementsSequence.find(g{2});
            ind = (indH - 1)*self.N.order + indN;
        end

        function s = computeSetN(self)
            s = replab.perm.Set.fromPermutationGroup(self.N.niceMorphism.target);
            s.sort;
        end

        function m = computeRegularMorphism(self)
            s = self.setN;
            o = double(self.N.order);
            N = self.N.niceMorphism.target; % use the permutation realization
            images = cell(1, N.nGenerators);
            for i = 1:N.nGenerators
                img = zeros(1, o);
                gen = N.generator(i);
                for j = 1:o
                    img(j) = s.find(gen(s.at(j)')');
                end
                images{i} = img;
            end
            m = self.N.morphismByImages(replab.S(o), 'images', images);
        end

        function a = computePhiMorphism(self)
            o = double(self.N.order);
            nGenH = self.H.nGenerators;
            phiImages = cell(1, nGenH);
            s = self.setN;
            for i = 1:nGenH
                phiImage = zeros(1, o);
                genH = self.H.generator(i);
                for j = 1:o
                    preimg_prm = s.at(j)';
                    preimg = self.N.niceMorphism.preimageElement(preimg_prm);
                    img = self.phi.leftAction(genH, preimg);
                    img_prm = self.N.niceMorphism.imageElement(img);
                    phiImage(j) = s.find(img_prm');
                end
                phiImages{i} = phiImage;
            end
            a = self.H.morphismByImages(replab.S(o), 'images', phiImages);
        end

    end

    methods (Access = protected) % Implementations

        % FiniteGroup

        function o = computeOrder(self)
            o = self.H.order * self.N.order;
        end

        function e = computeElementsSequence(self)
            e = replab.Sequence.lambda(self.order, ...
                                       @(ind) self.atFun(ind), ...
                                       @(g) self.findFun(g));
        end

        function s = computeSetProduct(self)
            TH = self.H.setProduct.sets;
            TN = self.N.setProduct.sets;
            idN = self.N.identity;
            idH = self.H.identity;
            TH1 = cellfun(@(t) cellfun(@(h) {h idN}, t, 'uniform', 0), TH, 'uniform', 0);
            TN1 = cellfun(@(t) cellfun(@(n) {idH n}, t, 'uniform', 0), TN, 'uniform', 0);
            s = replab.SetProduct(self, horzcat(TH1, TN1), true);
        end

    end

    methods % Implementations

        % Domain

        function s = sample(self)
            s = sample@replab.SemidirectProductGroup(self);
        end

        % FiniteSet

        function res = hasSameTypeAs(self, rhs)
            if isa(rhs.type, 'replab.SemidirectProduct') && isa(rhs.type, 'replab.FiniteGroup')
                if self.type.H == rhs.type.H && self.type.N == rhs.type.N
                    for i = 1:self.H.type.nGenerators
                        h = self.H.type.generator(i);
                        for j = 1:self.N.type.nGenerators;
                            n = self.N.type.generator(j);
                            if ~self.N.type.eqv(self.type.phi.leftAction(h, n), rhs.type.phi.leftAction(h, n))
                                res = false;
                                return
                            end
                        end
                    end
                    res = true;
                else
                    res = false;
                end
            else
                res = false;
            end
        end

        % FiniteGroup

        function G = withGeneratorNames(self, newNames)
            if isequal(self.generatorNames, newNames)
                G = self;
                return
            end
            G = replab.prods.SemidirectProduct_finite(self.phi, self.type, 'generatorNames', newNames);
        end

    end

    methods

        function s = setN(self)
        % Returns the sorted set of all permutations in the nice morphism image of the acted upon group
            s = self.cached('setN', @() self.computeSetN);
        end

        function p = niceImage(self, x)
            h = x{1};
            n = x{2};
            p1 = self.H.niceMorphism.imageElement(h);
            p2a = self.phiMorphism.imageElement(h);
            p2b = self.regularMorphism.imageElement(n);
            p2 = p2a(p2b);
            p = [p1 p2+length(p1)];
        end

        function m = regularMorphism(self)
            m = self.cached('regularMorphism', @() self.computeRegularMorphism);
        end

        function a = phiMorphism(self)
            a = self.cached('phiMorphism', @() self.computePhiMorphism);
        end

    end

end
