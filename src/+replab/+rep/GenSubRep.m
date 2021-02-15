classdef GenSubRep < replab.Obj
% Describes a generic subrepresentation of a finite dimensional representation
%
% The subrepresentation is generic in the sense it can be defined over a division ring. For example, a quaternionic-type subrepresentation
% is defined using quaternion injection and projection maps.
%
% This class is present mostly to cater for the different cases present when refining subrepresentations.
%
% Note that this class is unrelated to `.Rep` in the inheritance hierarchy, but it follows the naming conventions of `.SubRep` and
% implements a subset of its methods/properties.
%
% Main caveats, compared to `.Rep`, are:
%
% * No support for exact subrepresentations; injection and projection maps are of type 'double/sparse'
%
% * The matrices returned are possibly sparse

    properties (SetAccess = protected)
        parent % (`+replab.Rep`): Parent representation of dimension $D$
        group % (`+replab.CompactGroup`): Group being represented
        dimension % (integer): Representation dimension
        divisionRing % ('R', 'C', 'H'): Division ring over which the subrepresentation is defined
        isUnitary % (logical): Whether this subrepresentation is unitary
        injection % (double(D,d) or `.H`(D,d), may be sparse): Injection map
        projection % (double(d,D) or `.H`(d,D), may be sparse): Projection map
        mapsAreAdjoint % (logical): True if `.parent` is unitary and `.injection` is the conjugate transpose of `.projection`
    end

    methods (Static, Access = protected)

        function P = correctQuaternionPair(I, P)
        % Corrects a pair of injection/projection map after conversion from a quaternionic rep
        %
        % The pairs should diverge mostly up to a quaternion scalar, such that ``P * I`` is a block matrix
        % with mostly identical 4x4 blocks on the diagonal.
        %
        % This function will correct the projection map by inverting the average 4x4 block,
        % and additionally perform one Newton-Raphson iteration.
        %
        % Args:
        %   I (double(\*,\*)): Injection map
        %   P (double(\*,\*)): Projection map
        %
        % Returns:
        %   double(\*,\*): Corrected projection map
            replab.msg(2, 'Converting GenSubRep to SubRep, correcting scalar factor');
            replab.msg(2, '--------------------------------------------------------');
            d = size(I, 2);
            PI = P * I;
            replab.msg(2, 'Before correction, error = %6.2E', norm(PI - eye(d), 'fro'));
            % Compute the quaternion scalar block by averaging over the diagonal blocks
            block = PI(1:4,1:4);
            for i = 2:d/4
                range = (i-1)*4+(1:4);
                block = block + PI(range, range);
            end
            block = block/(d/4);
            % invert and correct
            P = kron(eye(d/4), inv(block)) * P;
            replab.msg(2, 'After block correction, error = %6.2E', norm(P*I - eye(d), 'fro'));
            % apply one Newton-Raphson iteration for increased precision
            P = 2*P - P*I*P;
            replab.msg(2, 'After Newton iteration, error = %6.2E', norm(P*I - eye(d), 'fro'));
        end

        function P = correctComplexPair(I, P)
        % Corrects a pair of injection/projection map after conversion from a complex rep
        %
        % The pairs should diverge mostly up to a complex scalar, such that ``P * I`` is a block matrix
        % with mostly identical 2x2 blocks on the diagonal.
        %
        % This function will correct the projection map by inverting the average 2x2 block,
        % and additionally perform one Newton-Raphson iteration.
        %
        % Args:
        %   I (double(\*,\*)): Injection map
        %   P (double(\*,\*)): Projection map
        %
        % Returns:
        %   double(\*,\*): Corrected projection map
            replab.msg(2, 'Converting GenSubRep to SubRep, correcting scalar factor');
            replab.msg(2, '--------------------------------------------------------');
            d = size(I, 2);
            PI = P * I;
            replab.msg(2, 'Before correction, error = %6.2E', norm(PI - eye(d), 'fro'));
            % Compute the complex scalar block by averaging over the diagonal blocks
            block = PI(1:2,1:2);
            for i = 2:d/2
                range = (i-1)*2+(1:2);
                block = block + PI(range, range);
            end
            block = block/(d/2);
            % invert and correct
            P = kron(eye(d/2), inv(block)) * P;
            replab.msg(2, 'After block correction, error = %6.2E', norm(P*I - eye(d), 'fro'));
            % apply one Newton-Raphson iteration for increased precision
            P = 2*P - P*I*P;
            replab.msg(2, 'After Newton iteration, error = %6.2E', norm(P*I - eye(d), 'fro'));
        end

    end

    methods (Access = protected)

        function sub = computeToSubRep(self)
            d = self.dimension;
            type = [self.divisionRing '->' self.parent.field]
            switch type
              case {'R->R', 'C->C'}
                if self.mapsAreAdjoint
                    sub = self.parent.subRep(self.injection, 'projection', self.projection, 'isUnitary', true);
                else
                    sub = self.parent.subRep(self.injection, 'projection', self.projection, 'isUnitary', false);
                end
              case 'C->R'
                A = real(self.injection);
                B = imag(self.injection);
                subI = [A, B]*sqrt(2);
                p = reshape(reshape(1:2*d, [d 2])', [1 d*2]);
                subI = subI(:, p);
                if self.mapsAreAdjoint
                    sub = self.parent.subRep(subI, 'projection', subI', 'divisionAlgebraName', 'complex', 'isUnitary', true);
                else
                    E = real(self.projection);
                    F = imag(self.projection);
                    subP = [E; -F]*sqrt(2);
                    subP = subP(p, :);
                    subP = replab.rep.GenSubRep.correctComplexPair(subI, subP);
                    sub = self.parent.subRep(subI, 'projection', subP, 'divisionAlgebraName', 'complex', 'isUnitary', false);
                end
              case 'H->R'
                A = self.injection.part1;
                B = self.injection.parti;
                C = self.injection.partj;
                D = self.injection.partk;
                d = self.dimension;
                p = reshape(reshape(1:4*d, [d 4])', [1 d*4]);
                subI = [A+B+C+D, A-B-C+D, A+B-C-D, A-B+C-D];
                subI = subI(:, p);
                if self.mapsAreAdjoint
                    sub = self.parent.subRep(subI, 'projection', subI', 'divisionAlgebraName', 'quaternion.rep', 'isUnitary', true);
                else
                    E = self.projection.part1;
                    F = self.projection.parti;
                    G = self.projection.partj;
                    H = self.projection.partk;
                    subP = [E-F-G-H; E+F+G-H; E-F+G+H; E+F-G+H];
                    subP = subP(p, :);
                    subP = replab.rep.GenSubRep.correctQuaternionPair(subI, subP);
                    sub = self.parent.subRep(subI, 'projection', subP, 'divisionAlgebraName', 'quaternion.rep', 'isUnitary', false);
                end
            end
        end

    end


    methods

        function self = GenSubRep(parent, divisionRing, isUnitary, injection, projection)
        % Constructs a generic subrepresentation
        %
        % Args:
        %   parent (`+replab.Rep`): Parent representation of dimension $D$
        %   divisionRing ('R', 'C', 'H'): Division ring over which the subrepresentation is defined
        %   isUnitary (logical): Whether this subrepresentation is unitary
        %   injection (double(D,d) or `.H`(D,d), may be sparse): Injection map
        %   projection (double(d,D) or `.H`(d,D), may be sparse): Projection map
            D = size(injection, 1);
            d = size(injection, 2);
            self.parent = parent;
            self.group = parent.group;
            self.dimension = d;
            self.divisionRing = divisionRing;
            self.isUnitary = isUnitary;
            self.injection = injection;
            self.projection = projection;
            self.mapsAreAdjoint = parent.knownUnitary && all(all(injection == projection'));
        end

        function sub = toSubRep(self)
        % Converts back the generic subrepresentation to a subrepresentation over the original field
        %
        % Returns:
        %   `+replab.SubRep`: Subrepresentation
            sub = self.cached('toSubRep', @() self.computeToSubRep);
        end

        function rho = image(self, g)
        % Computes the image of the given group element
        %
        % See `+replab.Rep.image`
            rho = self.projection * self.parent.image(g, 'double/sparse') * self.injection;
        end

        function rho = inverseImage(self, g)
        % Computes the image of the inverse of the given group element
        %
        % See `+replab.Rep.inverseImage`
            rho = self.projection * self.parent.inverseImage(g, 'double/sparse') * self.injection;
        end

        function rho = sample(self)
        % Returns a random image of this subrepresentation
        %
        % See `+replab.Rep.sample`
            rho = self.image(self.group.sample);
        end

        function mat = projector(self)
        % Returns a projector on the subspace of this subrepresentation
        %
        % Note: the projector can be complex/quaternion-valued.
        %
        % See `+replab.SubRep.projector`
            mat = self.injection * self.projection;
        end

        function gen1 = refine(self, varargin)
        % Refines a generic subrepresentation
        %
        % See `+replab.SubRep.refine`
            args = struct('numNonImproving', 20, 'largeScale', self.parent.dimension > 1000, 'nSamples', 5, 'nInnerIterations', 3, 'maxIterations', 1000, 'injectionBiortho', [], 'projectionBiortho', []);
            args = replab.util.populateStruct(args, varargin);
            if self.parent.knownUnitary
                if args.largeScale
                    gen1 = replab.rep.refine_unitary_largeScale(self, args.numNonImproving, args.nSamples, args.maxIterations, []);
                else
                    gen1 = replab.rep.refine_unitary_mediumScale(self, args.nInnerIterations, args.maxIterations, []);
                end
            else
                if args.largeScale
                    gen1 = replab.rep.refine_nonUnitaryLargeScale(self, args.numNonImproving, args.nSamples, args.maxIterations, [], []);
                else
                    gen1 = replab.rep.refine_nonUnitaryMediumScale(self, args.nInnerIterations, args.maxIterations, [], []);
                end
            end
        end

    end

    methods % Implementations

        function l = laws(self)
            l = replab.rep.GenSubRepLaws(self);
        end

    end

    methods (Static)

        function gen = fromSubRep(sub)
        % Creates a generic subrepresentation from a subrepresentation
        %
        % Exploits a division algebra if present (delegates to `.fromComplexTypeSubRep` or `.fromQuaternionTypeSubRep`).
        %
        % Guarantees that the roundtrip ``gen = replab.rep.GenSubRep.fromSubRep(sub); sub1 = gen.toSubRep`` preserves
        % the division algebra structure if present.
        %
        % Args:
        %   sub (`+replab.SubRep`): Subrepresentation
        %
        % Returns:
        %   `.GenSubRep`: Generic subrepresentation over the same field
            if sub.overR
                switch sub.divisionAlgebraName
                  case ''
                    gen = replab.rep.GenSubRep(sub.parent, 'R', sub.knownUnitary, sub.injection, sub.projection);
                  case 'complex'
                    gen = replab.rep.GenSubRep.fromComplexTypeSubRep(sub);
                  case 'quaternion.rep'
                    gen = replab.rep.GenSubRep.fromQuaternionTypeSubRep(sub);
                end
            else
                gen = replab.rep.GenSubRep(sub.parent, 'C', sub.knownUnitary, sub.injection, sub.projection);
            end
        end

        function gen = fromQuaternionTypeSubRep(sub)
        % Creates a generic subrepresentation from a real subrepresentation encoding a quaternion division algebra
        %
        % Args:
        %   sub (`+replab.SubRep`): Real subrepresentation
        %
        % Returns:
        %   `.GenSubRep`: Generic subrepresentation over the quaternions
            assert(sub.overR);
            i = replab.H.i;
            j = replab.H.j;
            k = replab.H.k;
            U = [1+i+j+k; 1-i+j-k; 1-i-j+k; 1+i-j-k]/4;
            d = sub.dimension;
            injection = sub.injection('double/sparse') * kron(eye(d/4), U);
            if sub.mapsAreAdjoint
                projection = injection';
                gen = replab.rep.GenSubRep(sub.parent, 'H', true, injection, projection);
            else
                projection = kron(eye(d/4), U)' * sub.projection('double/sparse');
                gen = replab.rep.GenSubRep(sub.parent, 'H', false, injection, projection);
            end
        end

        function gen = fromComplexTypeSubRep(sub)
        % Creates a generic subrepresentation from a real subrepresentation encoding a complex division algebra
        %
        % Args:
        %   sub (`+replab.SubRep`): Real subrepresentation
        %
        % Returns:
        %   `.GenSubRep`: Generic subrepresentation over the complex field
            assert(sub.overR);
            U = [1; 1i]/sqrt(2);
            d = sub.dimension;
            injection = sub.injection('double/sparse') * kron(eye(d/2), U);
            if sub.mapsAreAdjoint
                projection = injection';
                gen = replab.rep.GenSubRep(sub.parent, 'C', true, injection, projection);
            else
                projection = kron(eye(d/2), U)' * sub.projection('double/sparse');
                gen = replab.rep.GenSubRep(sub.parent, 'C', false, injection, projection);
            end
        end

    end

end
