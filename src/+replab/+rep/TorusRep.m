classdef TorusRep < replab.Rep
% Representation of a `+replab.TorusGroup` using diagonal matrices

    properties (SetAccess = protected)
        torusMap % (integer(\*,\*)): Morphism from the represented torus group to the torus used for the representation
    end

    methods

        function self = TorusRep(group, torusMap)
            assert(isa(group, 'replab.TorusGroup'));
            d = size(torusMap, 1);
            td = sum(all(torusMap == 0), 2);
            self@replab.Rep(group, 'C', d, 'isUnitary', true, 'trivialDimension', td);
            self.torusMap = torusMap;
        end

    end

    methods (Access = protected) % Implementations

        function e = computeErrorBound(self)
            e = sqrt(self.dimension)*1e-15; % conservative estimate
        end

        function rho = image_double_sparse(self, g)
            rho = replab.TorusGroup.torusRepImage(self.torusMap*g);
        end

    end

    methods % Implementations

        % Rep

        function b = isExact(self)
            b = false; % Note: redundant, as `.Rep.isExact` is already false
        end

        function p = invariantBlocks(self)
            p = replab.Partition.finest(self.dimension);
        end

        function b = hasTorusImage(self)
            b = true;
        end

        function [torusMap, torusInjection, torusProjection] = torusImage(self)
            d = size(self.torusMap, 1);
            torusMap = self.torusMap * self.group.injection;
            torusInjection = speye(d);
            torusProjection = speye(d);
        end

    end

    methods (Static)

        function [blocks1, blocks2] = matchTorusMaps(torusMap1, torusMap2)
        % Given two torus maps, matches the rows of both
        %
        % We have ``all(torusMap1(i1,:) == torusMap2(i2,:))`` for ``i1`` in ``blocks1(k)`` and ``i2`` in ``blocks2(k)``.
        %
        % Args:
        %   torusMap1 (integer(\*,n)): First torus map
        %   torusMap2 (integer(\*,n)): Second torus map
        %
        % Returns
        % -------
        %   blocks1: cell(1,\*) of integer(1,\*)
        %     Blocks of row indices in ``torusMap1``
        %   blocks2: cell(1,\*) of integer(1,\*)
        %     Blocks of row indices in ``torusMap2``
            b1 = find(all(torusMap1 == 0, 2));
            b2 = find(all(torusMap2 == 0, 2));
            if ~isempty(b1) && ~isempty(b2)
                blocks1 = {b1};
                blocks2 = {b2};
            else
                blocks1 = cell(1, 0);
                blocks2 = cell(1, 0);
            end
            r = 1;
            while r <= size(torusMap1, 1)
                row = torusMap1(r, :);
                if any(row ~= 0)
                    b1 = find(all(bsxfun(@eq, torusMap1, row), 2));
                    b2 = find(all(bsxfun(@eq, torusMap2, row), 2));
                    blocks1{1,end+1} = b1;
                    blocks2{1,end+1} = b2;
                    torusMap1(b1, :) = 0;
                    torusMap2(b2, :) = 0;
                end
                r = r + 1;
            end
        end

    end

end
