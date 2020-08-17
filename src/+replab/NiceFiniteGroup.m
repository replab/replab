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

        function G = niceGroup(self)
        % Returns the permutation group isomorphic to this finite group
        %
        % This is the image of `.niceMorphism`.
        %
        % Returns:
        %   `+replab.PermutationGroup`: The image of the isomorphism
            G = self.cached('niceGroup', @() self.computeNiceGroup);
        end

    end

    methods (Access = protected)

        function G = computeNiceGroup(self)
            gens = cellfun(@(g) self.niceImage(g), self.generators, 'uniform', 0);
            ds = length(self.niceImage(self.identity));
            G = replab.PermutationGroup(ds, gens, self.cachedOrEmpty('order'));
        end

        function order = computeOrder(self)
            order = self.niceGroup.order;
        end

        function m = computeNiceMorphism(self)
            m = replab.mrp.NiceFiniteGroupIsomorphism(self, self.niceGroup);
        end

        function dec = computeDecomposition(self)
            prmD = self.niceGroup.decomposition;
            T1 = cellfun(@(T) cellfun(@(t) self.niceMorphism.preimageElement(t), T, 'uniform', 0), prmD.T, 'uniform', 0);
            dec = replab.FiniteGroupDecomposition(self, T1);
        end

        function C = computeConjugacyClasses(self)
            C = cellfun(@(cl) replab.ConjugacyClass(self, self.niceMorphism.preimageElement(cl.representative), ...
                                                    self.niceMorphism.preimageGroup(cl.representativeCentralizer)), ...
                        self.niceGroup.conjugacyClasses, 'uniform', 0);
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
            niceGroup = self.niceGroup.subgroupWithGenerators(cellfun(@(g) self.niceImage(g), generators, 'uniform', 0), order);
            sub = self.niceSubgroup(generators, order, niceGroup);
        end

        function sub = niceSubgroup(self, generators, order, niceGroup)
            sub = replab.NiceFiniteSubgroup(self.type, generators, order);
            if nargin > 3 && ~isempty(niceGroup)
                sub.cache('niceGroup', niceGroup, '==');
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
            subgroup = self.niceMorphism.preimageGroup(B.subgroup);
            canRep = self.niceMorphism.preimageElement(B.representative);
            B = replab.LeftCoset(subgroup, canRep, self);
        end

        % Relation to other groups

        % Morphisms

        function m = morphismByImages(self, target, images)
            first = self.niceMorphism; % maps this to the perm group
            second = self.niceGroup.morphismByImages(target, images); % from the perm group to the images
            m = first.andThen(second);
        end

    end

end
