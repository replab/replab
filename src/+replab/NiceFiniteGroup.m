classdef NiceFiniteGroup < replab.FiniteGroup
% A nice finite group is a finite group equipped with an injective homomorphism into a permutation group
%
% Any class that subclasses `.NiceFiniteGroup` must implements a method `.niceImage` that returns a permutation
% corresponding to a group element.
%
% Then, an isomorphism is established between the present `.NiceFiniteGroup` and a permutation group; as
% permutation groups can be handled by efficient BSGS algorithms, the requested computations can be
% translated back and forth between this group and a permutation group.
%
% In particular, the decomposition of the finite group in a product of sets (`.FiniteGroupDecomposition`),
% the enumeration of elements using a `.IndexedFamily`, the construction of subgroups is all handled
% by permutation group algorithms.
%
% Note that the `.niceImage` method of this group should be valid for all elements of the parent as well,
% and that groups having the same type (as verified by `.hasSameTypeAs`) should return the same images
% under `.niceImage`.

    methods % Nice monomorphism support

        function g = nicePreimage(self, p)
        % Default implementation of chain
            g = self.niceMorphism.preimageElement(p);
        end

        function p = niceImage(self, g)
        % Returns a permutation image of the given group element
        %
        % Args:
        %   g (element): Group element to represent as a permutation
        %
        % Returns:
        %   permutation: Permutation image of ``g``
            error('Abstract');
        end

        function sub = niceSubgroup(self, generators, order, niceGroup)
        % Constructs a subgroup of this group with a small optional optimization
        %
        % This method is used by the other subgroup construction methods.
        %
        % Args:
        %   generators (cell(1,\*) of group elements): Subgroup generators
        %   order (vpi or ``[]``, optional): Subgroup order if known, default value ``[]``
        %   niceGroup (`.PermutationGroup` or ``[]``, optional): Permutation realization of the subgroup if known, default value ``[]``
            if nargin < 3 || isempty(order)
                order = [];
            end
            sub = replab.NiceFiniteSubgroup(self.type, generators, order);
            if nargin > 3 && ~isempty(niceGroup)
                sub.cache('niceGroup', niceGroup, '==');
            end
        end

    end

    methods (Access = protected) % Implementations

        function order = computeOrder(self)
            order = self.niceGroup.order;
        end

        function dec = computeDecomposition(self)
            prmD = self.niceGroup.decomposition;
            T1 = cellfun(@(T) cellfun(@(t) self.niceMorphism.preimageElement(t), T, 'uniform', 0), prmD.T, 'uniform', 0);
            dec = replab.FiniteGroupDecomposition(self, T1);
        end

        function C = computeConjugacyClasses(self)
            nc = self.niceGroup.conjugacyClasses;
            C = nc.imap(self.niceMorphism.inverse);
        end

        function res = computeIsCyclic(self)
            res = self.niceGroup.isCyclic;
        end

        function res = computeIsSimple(self)
            res = self.niceGroup.isSimple;
        end

        function E = computeElements(self)
            atFun = @(ind) self.niceMorphism.preimageElement(self.niceGroup.elements.at(ind));
            findFun = @(el) self.niceGroup.elements.find(self.niceImage(el));
            E = replab.IndexedFamily.lambda(self.order, atFun, findFun);
        end

        function sub = computeDerivedSubgroup(self)
            sub = self.niceMorphism.preimageGroup(self.niceGroup.derivedSubgroup);
        end

    end

    methods (Access = protected)

        function G = computeNiceGroup(self)
            gens = cellfun(@(g) self.niceImage(g), self.generators, 'uniform', 0);
            ds = length(self.niceImage(self.identity));
            G = replab.PermutationGroup(ds, gens, self.cachedOrEmpty('order'));
        end

        function m = computeNiceMorphism(self)
            m = replab.mrp.NiceFiniteGroupIsomorphism(self, self.niceGroup);
        end

        function A = computeDefaultAbstractGroup(self)
            A = self.niceGroup.abstractGroup;
        end

        function m = computeDefaultAbstractMorphism(self)
            m = self.niceMorphism.andThen(self.niceGroup.abstractMorphism);
        end

    end

    methods % Implementation

        % Domain

        function g = sample(self)
            g = self.niceMorphism.preimageElement(self.niceGroup.sample);
            % TODO: optimize
        end

        % FiniteSet

        % FiniteGroup

        % Group elements

        function b = contains(self, g)
            b = self.niceGroup.contains(self.niceImage(g));
        end

        function o = elementOrder(self, g)
            o = self.niceGroup.elementOrder(self.niceImage(g));
        end

        function l = factorizeLetters(self, element)
            l = self.niceGroup.factorizeLetters(self.niceImage(element));
        end

        % Construction of groups

        function res1 = closure(self, obj)
            if isa(obj, 'replab.FiniteGroup')
                % if one group contains the other
                if self.isSubgroupOf(obj)
                    res1 = obj;
                    return
                end
                if obj.isSubgroup(self)
                    res1 = self;
                    return
                end
                % otherwise do the computation
                res = self.niceGroup.closure(self.type.niceMorphism.imageGroup(obj));
            else
                % if the group already contains the element
                if self.contains(obj)
                    res1 = self;
                    return
                end
                % otherwise do the computation
                res = self.niceGroup.closure(self.type.niceImage(obj));
            end
            res1 = self.type.niceMorphism.preimageGroup(res);
        end

        function res1 = normalClosure(self, obj)
            if isa(obj, 'replab.FiniteGroup')
                res = self.niceGroup.normalClosure(self.type.niceMorphism.imageGroup(obj));
            else
                res = self.niceGroup.normalClosure(self.niceImage(obj));
            end
            res1 = self.type.niceMorphism.preimageGroup(res);
        end

        % Subgroups

        function sub = subgroupWithGenerators(self, generators, order)
            if nargin < 3
                order = [];
            end
            if length(generators) == self.nGenerators && all(arrayfun(@(i) self.eqv(self.generator(i), generators{i}), 1:length(generators)))
                sub = self;
            else
                niceGroup = self.niceGroup.subgroupWithGenerators(cellfun(@(g) self.niceImage(g), generators, 'uniform', 0), order);
                sub = self.niceSubgroup(generators, order, niceGroup);
            end
        end

        function sub1 = centralizer(self, obj)
            if isa(obj, 'replab.FiniteGroup')
                sub = self.niceGroup.centralizer(self.niceMorphism.imageGroup(obj));
            else
                sub = self.niceGroup.centralizer(self.niceImage(obj));
            end
            sub1 = self.niceMorphism.preimageGroup(sub);
        end

        function sub1 = intersection(self, rhs)
            sub1 = self.niceMorphism.preimageGroup(self.niceGroup.intersection(self.niceMorphism.imageGroup(rhs)));
        end

        % Cosets

        function B = findLeftConjugations(self, s, t, sCentralizer, tCentralizer)
            s = self.niceImage(s);
            t = self.niceImage(t);
            if nargin < 4 || isempty(sCentralizer)
                sCentralizer = [];
            else
                sCentralizer = self.niceMorphism.imageGroup(sCentralizer);
            end
            if nargin < 5 || isempty(tCentralizer)
                tCentralizer = [];
            else
                tCentralizer = self.niceMorphism.imageGroup(tCentralizer);
            end
            B = self.niceGroup.findLeftConjugations(s, t, sCentralizer, tCentralizer);
            group = self.niceMorphism.preimageGroup(B.group);
            canRep = self.niceMorphism.preimageElement(B.representative);
            B = replab.LeftCoset(group, canRep, self);
        end

        % Relation to other groups

    end

    methods (Access = protected)

        % Morphisms

        function m = morphismByImages_(self, target, preimages, images)
            first = self.niceMorphism; % maps this to the perm group
            preimagesNG = cellfun(@(g) self.niceMorphism.imageElement(g), preimages, 'uniform', 0);
            second = self.niceGroup.morphismByImages_(target, preimagesNG, images); % from the perm group to the images
            m = first.andThen(second);
        end

    end

end
