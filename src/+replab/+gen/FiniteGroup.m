classdef FiniteGroup < replab.FiniteGroup & replab.gen.FiniteSet
% A generic finite group

    methods

        function self = FiniteGroup(type, generators, varargin)
        % Constructs a nice finite group
        %
        % If ``nice`` is provided, then ``generatorNames``, ``order``, ``relators`` cannot be provided.
        % If ``nice`` is provided, then ``niceIsomorphism`` must be provided too.
        %
        % Args:
        %   type (`+replab.FiniteGroupType`): Finite group type
        %   generators (cell(1,\*) of group elements): Generators
        %
        % Keyword Args:
        %   generatorNames (cell(1,\*) of charstring, optional): Names of the generators
        %   nice (`+replab.FiniteGroup`): Nice group on which computations are performed
        %   niceIsomorphism (`+replab.+gen.NiceIsomorphism`): Isomorphism
        %   order (vpi or integer, optional): Group order
        %   relators (cell(1,\*) of charstring or integer(1,\*), optional): Relators
            args = struct('generatorNames', [], 'nice', [], 'niceIsomorphism', [], 'order', [], 'relators', []);
            args = replab.util.populateStruct(args, varargin);
            if isempty(args.nice)
                if isempty(args.niceIsomorphism)
                    niceIsomorphism = type.constructNiceIsomorphism(generators);
                else
                    niceIsomorphism = args.niceIsomorphism;
                end
                niceArgs = cell(1, 0);
                if ~isempty(args.generatorNames)
                    niceArgs{1,end+1} = 'generatorNames';
                    niceArgs{1,end+1} = args.generatorNames;
                end
                if ~isempty(args.order)
                    niceArgs{1,end+1} = 'order';
                    niceArgs{1,end+1} = args.order;
                end
                if ~isempty(args.relators)
                    niceArgs{1,end+1} = 'relators';
                    niceArgs{1,end+1} = args.relators;
                end
                niceGenerators = cellfun(@(g) niceIsomorphism.imageElement(g), generators, 'uniform', 0);
                nice = niceIsomorphism.target.type.groupWithGenerators(niceGenerators, niceArgs{:});
            else
                assert(~isempty(args.niceIsomorphism));
                assert(isempty(args.generatorNames));
                assert(isempty(args.order));
                assert(isempty(args.relators));
                nice = args.nice;
                niceIsomorphism = args.niceIsomorphism;
            end
            self@replab.gen.FiniteSet(type, nice, niceIsomorphism);
            self.identity = type.identity;
            self.representative_ = type.identity;
            self.generators = generators;
        end

    end

    methods (Access = protected) % Implementations

        % Morphisms

        function m = morphismByImages_(self, target, preimages, images, imageElementFun)
            first = self.niceIsomorphism; % maps this to the perm group
            preimagesNG = cellfun(@(g) first.imageElement(g), preimages, 'uniform', 0);
            second = self.nice.morphismByImages(target, 'preimages', preimagesNG, 'images', images, 'nChecks', 0); % from the perm group to the images
            if isa(second, 'replab.FiniteMorphism')
                m = replab.mrp.FiniteComposition(second, first, imageElementFun);
            else
                m = replab.mrp.Composition(second, first, imageElementFun);
            end
        end

    end

    methods % Implementation

        % Str

        function s = headerStr(self)
            if self.knownOrder
                s = sprintf('Finite group of order %s', strtrim(num2str(self.order)));
            else
                s = 'Finite group';
            end
        end

        function names = hiddenFields(self)
            names1 = hiddenFields@replab.gen.FiniteSet(self);
            names2 = hiddenFields@replab.FiniteGroup(self);
            names = union(names1, names2);
        end

        % Domain

        function b = eqv(self, x, y)
            b = self.type.eqv(x, y);
        end

        function g = sample(self)
            g = self.niceIsomorphism.preimageElement(self.nice.sample);
        end

        % Monoid

        function z = compose(self, x, y)
            z = self.type.compose(x, y);
        end

        % Group

        function y = inverse(self, x)
            y = self.type.inverse(x);
        end

        % FiniteGroup

        function a = abelianInvariants(self)
            a = self.nice.abelianInvariants;
        end

        function c1 = centralizer(self, obj)
            if isa(obj, 'replab.FiniteGroup')
                c = self.nice.centralizer(self.niceIsomorphism.imageGroup(obj));
            else
                c = self.nice.centralizer(self.niceIsomorphism.imageElement(obj));
            end
            c1 = self.niceIsomorphism.preimageGroup(c);
        end

% $$$         function res = closure()
% $$$             TODO
% $$$         end

% $$$         function res1 = closure(self, obj)
% $$$             if isa(obj, 'replab.FiniteGroup')
% $$$                 % if one group contains the other
% $$$                 if self.isSubgroupOf(obj)
% $$$                     res1 = obj;
% $$$                     return
% $$$                 end
% $$$                 if obj.isSubgroup(self)
% $$$                     res1 = self;
% $$$                     return
% $$$                 end
% $$$                 % otherwise do the computation
% $$$                 % TODO: enlarge isomorphism if necessary
% $$$                 res = self.nice.closure(self.niceIsomorphism.imageGroup(obj));
% $$$             else
% $$$                 % if the group already contains the element
% $$$                 if self.contains(obj)
% $$$                     res1 = self;
% $$$                     return
% $$$                 end
% $$$                 % otherwise do the computation
% $$$                 res = self.niceGroup.closure(self.type.niceImage(obj));
% $$$             end
% $$$             res1 = self.type.niceMorphism.preimageGroup(res);
% $$$         end


% $$$         function c = conjugacyClass(self, el, varargin)
% $$$             TODO
% $$$         end

% $$$         function C = conjugacyClass(self, varargin)
% $$$         % TODO
% $$$     end

        function sub = derivedSubgroup(self)
            sub = self.niceIsomorphism.preimageGroup(self.nice.derivedSubgroup);
        end

        function o = elementOrder(self, g)
            o = self.nice.elementOrder(self.niceIsomorphism.imageElement(g));
        end

        function e = exponent(self)
            e = self.nice.exponent;
        end

        function l = factorizeFlat(self, element)
            if isa(element, 'replab.FiniteSet')
                l = self.nice.factorizeFlat(element.imap(self.niceIsomorphism));
            else
                l = self.nice.factorizeFlat(self.niceIsomorphism.imageElement(element));
            end
        end

        function C = findLeftConjugations(self, s, t, varargin)
            mapGroup = @(G) self.niceIsomorphism.imageGroup(G);
            args1 = replab.kv.map(varargin, struct('sCentralizer', mapGroup, 'tCentralizer', mapGroup));
            s1 = self.niceIsomorphism.imageElement(s);
            t1 = self.niceIsomorphism.imageElement(t);
            C1 = self.nice.findLeftConjugations(s1, t1, args{:});
            C = C1.imap(self.niceIsomorphism.inverse);
        end

        function names = generatorNames(self)
            names = self.nice.generatorNames;
        end

        % TODO: rhs is not necessarily contained in this group! so need to construct
        % a larger closure group first in that case
% $$$         function sub1 = intersection(self, rhs)
% $$$             sub1 = self.niceMorphism.preimageGroup(self.niceGroup.intersection(self.niceMorphism.imageGroup(rhs)));
% $$$         end

        function res = isCyclic(self)
            res = self.nice.isCyclic;
        end

        function res = isSimple(self)
            res = self.nice.isSimple;
        end

        function res = knownOrder(self)
            res = self.nice.knownOrder;
        end

        % TODO: enlarge the isomorphism if necessary
% $$$         function res1 = normalClosure(self, obj)
% $$$             if isa(obj, 'replab.FiniteGroup')
% $$$                 res = self.niceGroup.normalClosure(self.type.niceMorphism.imageGroup(obj));
% $$$             else
% $$$                 res = self.niceGroup.normalClosure(self.type.niceImage(obj));
% $$$             end
% $$$             res1 = self.type.niceMorphism.preimageGroup(res);
% $$$         end

        function o = order(self)
            o = self.nice.order;
        end



        function R = relatorsFlat(self)
            R = self.nice.relatorsFlat;
        end

% $$$         function rep = regularRep(self)
% $$$             rep = self.niceMorphism.andThen(self.niceGroup.regularRep);
% $$$         end


% $$$         function res = withGeneratorNames(self, newNames)
% $$$         % TODO
% $$$         end

        % FiniteGroup/Cosets

        % doubleCoset / doubleCosets / leftCoset / leftCosets / normalCoset / normalCosets
        % rightCoset / rightCosets

        % FiniteGroup / Morphism

% $$$         function A = abstractGroup(self, generatorNames)
% $$$             if nargin < 2 || isempty(generatorNames)
% $$$                 generatorNames = self.generatorNames;
% $$$             end
% $$$             A = self.genericGroup.abstractGroup(generatorNames);
% $$$         end
% $$$
% $$$         function m = abstractMorphism(self, generatorNames)
% $$$             if nargin < 2 || isempty(generatorNames)
% $$$                 generatorNames = self.generatorNames;
% $$$             end
% $$$             m = self.genericIsomorphism.andThen(self.genericGroup.abstractMorphism(generatorNames));
% $$$         end

        % permutationIsomorphism

        % Subgroups

% $$$

        % Cosets

    end

end
