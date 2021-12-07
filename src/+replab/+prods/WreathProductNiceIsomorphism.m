classdef WreathProductNiceIsomorphism < replab.gen.NiceIsomorphism

    properties (SetAccess = protected)
        H % (`+replab.PermutationGroup`): Permutation group
        A % (`+replab.FiniteGroup`): Factor group
        isomorphism % (`+replab.FiniteIsomorphism`): Isomorphism from `.A` to a permutation group
        n % (integer): Number of copies of the factor group
    end

    methods

        function self = WreathProductNiceIsomorphism(type, H, A)
            assert(isa(A, 'replab.FiniteGroup'));
            n = type.n;
            self.H = H;
            self.A = A;
            self.n = n;
            isomorphism = A.orderPreservingPermutationIsomorphism;
            self.isomorphism = isomorphism;
            d = isomorphism.target.domainSize;
            % source/target construction
            order = self.H.order * A.order^n;
            sourceGenerators = cellfun(@(h) {h, type.identity{2}}, H.generators, 'uniform', 0);
            generatorNames = H.generatorNames;
            for i = 1:n
                for j = 1:A.nGenerators
                    g = type.identity;
                    g{2}{i} = A.generator(j);
                    sourceGenerators{1,end+1} = g;
                    generatorNames{1,end+1} = sprintf('%s%d', A.generatorNames{j}, i);
                end
            end
            targetGenerators = cellfun(@(g) self.imageElement(g), sourceGenerators, 'uniform', 0);
            if length(unique(generatorNames)) ~= length(generatorNames) % names are not unique
                generatorNames = [];
            end
            nice = replab.PermutationGroup(n*d, targetGenerators, 'order', order, 'generatorNames', generatorNames);
            self.source = replab.prods.WreathProductGroup_finite(H, A, type, sourceGenerators, nice, self);
            self.target = nice;
        end

    end

    methods % Implementations

        % Morphism

        function t = imageElement(self, s)
            n = self.n;
            h = s{1};
            base = s{2};
            iso = self.isomorphism;
            im = iso.imageElement(base{1});
            d = length(im);
            basePerm = im;
            shift = d;
            for i = 2:n
                im = iso.imageElement(base{i}) + shift;
                basePerm = [basePerm im];
                shift = shift + d;
            end
            ip = reshape(1:n*d, [d n]);
            ip = ip(:,h);
            ip = ip(:)';
            t = ip(basePerm);
        end

        % Isomorphism

        function s = preimageElement(self, t)
            n = self.n;
            d = self.isomorphism.target.domainSize;
            t = reshape(t, [d n]);
            h = (min(t, [], 1) - 1)/d + 1;
            ip = reshape(1:n*d, [d n]);
            ip(:,h) = ip;
            ip = ip(:)';
            t = ip(t);
            base = cell(1, n);
            for i = 1:n
                base{i} = self.isomorphism.preimageElement(t(:,i)'-(i-1)*d);
            end
            s = {h, base};
        end

        % NiceIsomorphism

        function l = sourceContains(self, s)
            l = self.H.contains(s{1}) && all(cellfun(@(a) self.A.contains(a), s{2}));
        end

    end

end
