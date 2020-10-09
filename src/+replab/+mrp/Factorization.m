classdef Factorization < replab.Obj
% Computes the factorization of a permutation group elements in its generators

    properties (SetAccess = protected)
        group % (`+replab.PermutationGroup`): Group which elements are factorized
        generators % (cell(1,\*) of elements of `.group`): Generators on which to factorize group elements
        useInverses % (logical): Whether the factorization can include generator inverses
    end

    methods

        function letters = factorize(self, g)
        % Returns the letters that compose the word with the given image
        %
        % Args:
        %   g (permutation in `.group`): Permutation to factorize
        %
        % Returns:
        %   integer(1,\*): Corresponding word expressed in letters
            error('Abstract');
        end

    end

    methods (Static)

        function f = make(group, generators, useInverses)
        % Constructs a group factorization object
        %
        % Args:
        %   group (`+replab.PermutationGroup`): Group to decompose elements of
        %   generators (cell(1,\*) of elements of ``group``): Group generators (default: ``group.generators``)
        %   useInverses (logical, optional): Whether to use inverses in the decomposition (default: true)
        %
        % Returns:
        %   `+replab.+mrp.Factorization`: The factorization object that can compute preimages in words of the generators
            if nargin < 2 || isempty(generators)
                generators = group.generators;
            end
            if nargin < 3 || isempty(useInverses)
                useInverses = true;
            end
            if group.isTrivial % treat trivial case
                f = replab.mrp.FactorizationTrivial(group, useInverses);
            elseif group.order <= replab.globals.factorizationOrderCutoff
                f = replab.mrp.FactorizationEnumeration.make(group, generators, useInverses);
            else
                f = replab.mrp.FactorizationChain(group, generators, useInverses);
            end
        end

    end

end
