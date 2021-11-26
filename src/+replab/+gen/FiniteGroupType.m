classdef FiniteGroupType < replab.FiniteGroupType
% Finite group type that delegates computations through an isomorphism
%
% Every object created in this group type has an isomorphic counterpart, its "nice"
% counterpart, stored in `+replab.+gen.FiniteSet.nice`.
%
% The classes defined in the `+replab.+gen` package solve two challenges:
%
% * Construct an isomorphism for any generating set of type elements.
%   Different strategies are used depending on the type of group considered.
%
% * Implement the different finite objects in RepLAB through this isomorphism scheme.
%   The implementations in `+replab.+gen.FiniteSet`, `+replab.+gen.LeftCoset`, ...
%   are straightforward.
%
% This class has an important subclass, `+replab.+gen.StaticFiniteGroupType`, used for
% group types where the same isomorphism can be reused for all generating sets of elements.
% For example, all signed permutation groups acting on the same domain can be mapped to
% the permutation group (which acts on a domain of doubled size).
%
% Other group types, such as the one for matrix groups with cyclotomic coefficients
% (see `+replab.+matrix.FiniteGroupType` ), construct tailor-made isomorphisms for
% each finite matrix group.

    methods (Access = protected)

        function G = groupFromNiceImage_(self, generators, nice, niceIsomorphism)
        % Creates a generic group of this type
        %
        % See `.constructGroup`, here the ``generators`` argument has been reconstructed if necessary
            G = replab.gen.FiniteGroup(self, generators, 'nice', nice, 'niceIsomorphism', niceIsomorphism);
        end

    end

    methods

        function iso = niceIsomorphism(self, elements)
        % Constructs a generic isomorphism for the given elements
        %
        % It is not guaranteed that the type of the isomorphism target is always the same.
        %
        % Args:
        %   elements (cell(1,\*) of group type elements): Elements (cannot be groups!)
        %
        % Returns:
        %   `.NiceIsomorphism`: An order preserving isomorphism whose source contains all given elements
            error('Abstract');
        end

        function [iso, varargout] = niceImages(self, varargin)
        % Finds images of the given groups/elements through a common isomorphism
        %
        % The method returns the nice isomorphism, and the image of all arguments through that isomorphism.
        %
        % The arguments can be either elements of this type or finite groups of this type.
        %
        % Returns
        % -------
        % iso: `.NiceIsomorphism`
        %   Order preserving isomorphism whose source contains all given elements
        % varargin:
        %   Images of the arguments
            n = length(varargin);
            areGroups = cellfun(@(x) isa(x, 'replab.FiniteGroup'), varargin);
            groups = varargin(areGroups);
            iso = [];
            % find the group with maximal order, if there is one, and reuse its nice isomorphism if it is compatible
            if ~isempty(groups)
                orders = cellfun(@(x) x.order, groups);
                ind = find(orders == max(orders), 1);
                iso = groups{ind}.niceIsomorphism;
                for i = 1:n
                    arg = varargin{i};
                    if isa(arg, 'replab.FiniteGroup')
                        if ~arg.compatibleWithNiceIsomorphism(iso)
                            iso = [];
                            break
                        end
                    else % group element
                        if ~iso.sourceContains(arg)
                            iso = [];
                            break
                        end
                    end
                end
            end
            if isempty(iso)
                % no preexisting isomorphism compatible with all elements, construct one
                elements = cell(1, 0);
                for i = 1:n
                    arg = varargin{i};
                    if isa(arg, 'replab.FiniteGroup')
                        elements = horzcat(elements, arg.generators);
                    else
                        elements{1,end+1} = arg;
                    end
                end
                iso = self.niceIsomorphism(elements);
            end
            % now we have a nice isomorphism, let's map the elements
            varargout = cell(1, n);
            for i = 1:n
                arg = varargin{i};
                res = [];
                if isa(arg, 'replab.FiniteGroup')
                    if arg.compatibleWithNiceIsomorphism(iso)
                        res = arg.nice;
                    else
                        res = iso.imageGroup(arg);
                    end
                else
                    res = iso.imageElement(arg);
                end
                varargout{i} = res;
            end
        end

        function G = groupFromNiceImage(self, nice, niceIsomorphism, varargin)
        % Creates a generic group of this type from its isomorphism image
        %
        % Args:
        %   nice (`+replab.FiniteGroup`): Subgroup of the isomorphism target whose generators are in 1-to-1 correspondance with ``generators``
        %   niceIsomorphism (`+replab.+gen.NiceIsomorphism`): Isomorphism whose source contains all the generators
        %
        % Keyword Args:
        %   generators (cell(1,\*) of elements of this type): Group generators
        %
        % Returns:
        %   `.FiniteGroup`: Constructed finite group
            args = struct('generators', []);
            args = replab.util.populateStruct(args, varargin);
            if isempty(args.generators)
                generators = cellfun(@(g) niceIsomorphism.preimageElement(g), nice.generators, 'uniform', 0);
            else
                generators = args.generators;
            end
            G = self.groupFromNiceImage_(generators, nice, niceIsomorphism);
        end

    end

    methods % Implementations

        function G = groupWithGenerators(self, generators, varargin)
            assert(all(ismember(varargin(1:2:end), {'generatorNames', 'order', 'relators'})));
            niceIso = self.niceIsomorphism(generators);
            target = niceIso.target;
            targetGenerators = cellfun(@(g) niceIso.imageElement(g), generators, 'uniform', 0);
            nice = target.subgroup(targetGenerators, varargin{:});
            G = self.groupFromNiceImage_(generators, nice, niceIso);
        end

    end

end
