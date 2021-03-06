classdef SequentialRep < replab.Rep
% A representation which consists of the sequential action of representations
%
% We write ``image(g) = reps{1}.image(g) * ... * reps{n}.image(g)``. The representations applied in sequence should be
% constructed in a way the product is a representation. For example, the construction works if the representations in
% the sequence commute.

    properties (SetAccess = protected)
        reps % (cell(1,\*) of `+replab.Rep`): Representations
    end

    methods

        function self = SequentialRep(group, field, dimension, reps)
        % Constructs a representation from a sequence of representations
        %
        % All the representations in the sequence should be defined on the same group, and on the same field.
        %
        % Args:
        %   group (`+replab.CompactGroup`): Common group
        %   field ({'R', 'C'}): Real or complex field
        %   dimension (integer): Representation dimension
        %   reps (cell(1,\*) of `+replab.Rep`): Representations
            replab.rep.assertCompatibleFactors(group, field, reps);
            assert(all(cellfun(@(r) r.dimension, reps) == dimension));
            isUnitary = cellfun(@(r) r.isUnitary, reps);
            self@replab.Rep(group, field, dimension, 'isUnitary', all(isUnitary));
            self.reps = reps;
        end

    end

    methods

        function n = nReps(self)
        % Returns the number of representations in the sequence
        %
        % Returns:
        %   integer: Number of representations
            n = length(self.reps);
        end

        function r = rep(self, i)
        % Returns the i-th representation in the sequence
        %
        % Args:
        %   i (integer): Index of the representation
        %
        % Returns:
        %   `+replab.Rep`: Representation
            r = self.reps{i};
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function c = decomposeTerm(self)
            c = self.reps;
        end

        function r = composeTerm(self, newFactors)
            r = replab.rep.SequentialRep(self.group, self.field, self.dimension, newFactors);
        end

        function rho = image_exact(self, g)
            if self.nReps == 0
                rho = replab.cyclotomic.eye(self.dimension);
            else
                rho = self.rep(1).image(g, 'exact');
                for i = 2:self.nReps
                    rho = rho * self.rep(i).image(g, 'exact');
                end
            end
        end

        function rho = image_double_sparse(self, g)
            if self.nReps == 0
                rho = speye(self.dimension);
            else
                rho = self.rep(1).image(g, 'double/sparse');
                for i = 2:self.nReps
                    rho = rho * self.rep(i).image(g, 'double/sparse');
                end
            end
        end

        function e = computeErrorBound(self)
            E = cellfun(@(r) r.errorBound, self.reps);
            C = cellfun(@(r) r.conditionNumberEstimate, self.reps);
            % we compute
            % e = e1*c2*...*cn + c1*e2*...*cn + ...
            e = 0;
            for i = 1:self.nReps
                e = e + E(i)*prod(C(1:i-1))*prod(C(i+1:end));
            end
        end

        function c = computeConditionNumberEstimate(self)
            c = prod(cellfun(@(r) r.conditionNumberEstimate, self.reps));
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function M = matrixRowAction_double_sparse(self, g, M)
            for i = self.nReps:-1:1
                M = self.rep(i).matrixRowAction(g, M, 'double/sparse');
            end
        end

        function M = matrixColAction_double_sparse(self, g, M)
            for i = 1:self.nReps
                M = self.rep(i).matrixColAction(g, M, 'double/sparse');
            end
        end

    end

    methods % Implementations

        % Str

        function names = hiddenFields(self)
            names = hiddenFields@replab.Rep(self);
            names{1, end+1} = 'reps';
        end

        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Rep(self);
            for i = 1:self.nReps
                names{1, end+1} = sprintf('rep(%d)', i);
                values{1, end+1} = self.rep(i);
            end
        end

        function s = headerStr(self)
            s = headerStr@replab.Rep(self); % logic in parent class
        end

        % Rep

        function b = isExact(self)
            b = all(cellfun(@(r) r.isExact, self.reps));
        end

        function p = invariantBlocks(self)
            if self.nReps == 0
                p = replab.Partition.trivial(self.dimension);
            else
                p = self.rep(1).invariantBlocks;
                for i = 2:self.nReps
                    p = p.join(self.rep(i).invariantBlocks);
                end
            end
        end

        function b = hasTorusImage(self)
            if any(cellfun(@(r) ~r.hasTorusImage, self.reps))
                b = false;
                return
            end
            nonTrivial = find(cellfun(@(r) any(any(r.torusImage ~= 0)), self.reps));
            if length(nonTrivial) <= 1
                b = true;
            else
                [~, ti1, ~] = self.rep(nonTrivial(1)).torusImage;
                for i = nonTrivial(2:end)
                    [~, ti, ~] = self.rep(i).torusImage;
                    if ~all(all(ti == ti1))
                        b = false;
                        return
                    end
                end
                b = true;
            end
        end

        function [torusMap, torusInjection, torusProjection] = torusImage(self)
            nonTrivial = find(cellfun(@(r) any(any(r.torusImage ~= 0)), self.reps));
            r = self.group.reconstruction.source.n; % torus rank
            if isempty(nonTrivial)
                torusMap = zeros(self.dimension, r);
                torusInjection = speye(self.dimension);
                torusProjection = speye(self.dimension);
            else
                [torusMap, torusInjection, torusProjection] = self.rep(nonTrivial(1)).torusImage;
                for i = nonTrivial(2:end)
                    torusMap = torusMap + self.rep(i).torusImage;
                end
            end
        end

    end

end
