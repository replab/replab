classdef FiniteFPGroup < replab.NiceFiniteGroup
% Describes a finite finitely presented group
%
% The group is the quotient of a free group by the normal closure of a set of relators. The relators are
% written using words in the free group.

    properties (SetAccess = protected)
        freeGroup % (`.FreeGroup`): Free group over which this group is defined
        relators % (cell(1,\*) of `.FreeGroupWord`): Relators
    end

    methods

        function self = FiniteFPGroup(freeGroup, relators)
        % Creates a finite finitely presented group from a free group and relators
        %
        % Example:
        %   >>> [F, a, x] = replab.FreeGroup.of('a', 'x');
        %   >>> replab.FiniteFPGroup(F, {a*a, x*x, a*x*inv(a)*inv(x)})
            assert(all(cellfun(@(r) r.group == freeGroup, relators)));
            self.type = self;
            self.freeGroup = freeGroup;
            self.relators = relators;
            self.identity = replab.FiniteFPGroupElement(self, freeGroup.identity);
            self.generators = cellfun(@(l) replab.FiniteFPGroupElement(self, l), freeGroup.generators, 'uniform', 0);
        end

    end

    methods (Static)

        function [fpg, varargout] = of(freeGroup, varargin)
        % Creates a finitely presented group from a description string
        %
        % Returns the finite group generators as additional output arguments.
        %
        % Example:
        %   >>> [F, x_] = replab.FreeGroup.of('x');
        %   >>> [G, x] = replab.FiniteFPGroup.of(F, x_*x_);
        %
        % Args:
        %   str (charstring): Single-line description string
        %
        % Returns:
        %   `+replab.FiniteFPGroup`: The parsed finitely presented group
            fpg = replab.FiniteFPGroup(freeGroup, varargin);
            if nargout > 1
                for i = 1:fpg.nGenerators
                    varargout{i} = fpg.generator(i);
                end
            end
        end

        function [fpg, varargout] = parsePresentation(str)
        % Creates a finitely presented group from a description string
        %
        % Returns the finite group generators as additional output arguments.
        %
        % Example:
        %   >>> [G, x] = replab.FiniteFPGroup.parsePresentation('< x | x^3 = 1 >');
        %
        % Args:
        %   str (charstring): Single-line description string
        %
        % Returns:
        %   `+replab.FiniteFPGroup`: The parsed finitely presented group
            P = replab.fp.Parser;
            [ok, names, relatorLetters] = P.parsePresentation(str);
            assert(ok, 'Error in given presentation string');
            F = replab.FreeGroup(names);
            % reduce the relators
            relators = cellfun(@(r) F.word(r), relatorLetters, 'uniform', 0);
            mask = cellfun(@(r) r.isEmpty, relators);
            % remove empty relators
            relators = relators(~mask);
            fpg = replab.FiniteFPGroup(F, relators);
            if nargout > 1
                for i = 1:length(fpg.nGenerators)
                    varargout{i} = fpg.generator(i);
                end
            end
        end

    end

    methods

        function nR = nRelators(self)
        % Returns the number of relators describing this group
        %
        % Returns:
        %   integer: Number of relators
            nR = length(self.relators);
        end

        function r = relator(self, i)
        % Returns the ``i``-th relator in the presentation of this group
        %
        % Args:
        %   i (integer): Index of relator
        % Returns:
        %   `+replab.FreeGroupWord`: ``i``-th relator
            r = self.relators{i};
        end

        function s = presentation(self)
        % Returns a readable text presentation of this group
        %
        % It can be parsed using `+replab.FiniteFPGroup.parsePresentation`.
        %
        % Returns:
        %   charstring: Presentation as string
            gens = strjoin(cellfun(@(e) e.toString, self.freeGroup.generators, 'uniform', 0), ', ');
            rels = strjoin(cellfun(@(r) r.toString, self.relators, 'uniform', 0), ' = ');
            s = ['< ' gens ' | ' rels ' = 1 >'];
        end

        function setPermutationImages(self, permutations)
        % Sets the permutation realization of this group
        %
        % Enables to skip the Todd-Coxeter procedure when operating over the group.
            self.cache('niceGroup', replab.PermutationGroup.of(permutations{:}), 'ignore');
            for i = 1:length(self.relators)
                assert(self.niceGroup.isIdentity(self.niceImage(replab.FiniteFPGroupElement(self, self.relators{i}))));
            end
        end

        function l = imagesDefineMorphism(self, target, generatorImages)
        % Checks whether the given images satisfy the relations of the presentation of this group
        %
        % If it returns true, it means those images describe a valid homomorphism from this `.FiniteFPGroup`
        % to the given target group.
        %
        % Args:
        %   target (`+replab.Group`): Target group
        %   generatorImages (cell(1,\*) of elements of ``target``): Images of the generators of this group
        %
        % Returns:
        %   logical: True if the images verify the presentation
            nR = length(self.relators);
            for i = 1:nR
                r = self.relators{i};
                if ~target.isIdentity(r.computeImage(target, generatorImages))
                    l = false;
                    return
                end
            end
            l = true;
        end

    end

    methods % Implementations

        function res = ne(self, rhs)
            res = ~(self == rhs);
        end

        function res = eq(self, rhs)
            if self.freeGroup ~= rhs.freeGroup
                res = false;
                return
            end
            % TODO: optimize
            allRelators = horzcat(self.relators, rhs.relators);
            fpg = self.freeGroup / allRelators;
            res = (self.order == fpg.order);
        end

        % Str

        function s = headerStr(self)
            s = self.presentation;
        end

        % Obj

        function l = laws(self)
            l = replab.FiniteFPGroupLaws(self);
        end

        % Domain

        function res = eqv(self, x, y)
            res = self.niceGroup.eqv(self.niceImage(x), self.niceImage(y));
        end

        % Monoid

        function z = compose(self, x, y)
            z = replab.FiniteFPGroupElement(self, x.representative * y.representative);
        end

        % Group

        function z = inverse(self, x)
            z = replab.FiniteFPGroupElement(self, inv(x.representative));
        end

        % NiceFiniteGroup

        function p = niceImage(self, x)
            p = x.representative.computeImage(self.niceGroup, self.niceGroup.generators);
        end

        % TODO: optimize morphismByImages imageElement

    end


end
