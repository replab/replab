classdef GenSubRep < replab.Obj
% Describes a generic subrepresentation of a finite dimensional representation
%
% The subrepresentation is generic in the sense it can be defined over a division ring. For example, a quaternionic-type subrepresentation
% is defined using quaternion injection and projection maps.
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
        parent % (`.Rep`): Parent representation of dimension $D$
        group % (`.CompactGroup`): Group being represented
        dimension % (integer): Representation dimension
        divisionRing % ('R', 'C', 'H'): Division ring over which the subrepresentation is defined
        isUnitary % (logical): Whether this subrepresentation is unitary
        injection % (double(D,d) or `.H`(D,d), may be sparse): Injection map
        projection % (double(d,D) or `.H`(d,D), may be sparse): Projection map
        mapsAreAdjoint % (logical): True if `.parent` is unitary and `.injection` is the conjugate transpose of `.projection`
    end

    methods (Static, Access = protected)

        function [I, P] = correctQuaternionPair(I, P)
            replab.msg(2, 'Converting GenSubRep to SubRep, correcting scalar factor');
            replab.msg(2, '--------------------------------------------------------');
            d = size(I, 2);
            PI = P * I;
            replab.msg(2, 'Before correction, error = %6.2E', norm(PI - eye(d), 'fro'));
            block = PI(1:4,1:4);
            for i = 2:d/4
                range = (i-1)*4+(1:4);
                block = block + PI(range, range);
            end
            block = block/(d/4);
            P = kron(eye(d/4), inv(block)) * P;
            replab.msg(2, 'After block correction, error = %6.2E', norm(P*I - eye(d), 'fro'));
            % one Newton-Raphson iteration
            P = 2*P - P*I*P;
            replab.msg(2, 'After Newton iteration, error = %6.2E', norm(P*I - eye(d), 'fro'));
        end

        function [I, P] = correctComplexPair(I, P)
            replab.msg(2, 'Converting GenSubRep to SubRep, correcting scalar factor');
            replab.msg(2, '--------------------------------------------------------');
            d = size(I, 2);
            PI = P * I;
            replab.msg(2, 'Before correction, error = %6.2E', norm(PI - eye(d), 'fro'));
            block = PI(1:2,1:2);
            for i = 2:d/2
                range = (i-1)*2+(1:2);
                block = block + PI(range, range);
            end
            block = block/(d/2);
            P = kron(eye(d/2), inv(block)) * P;
            replab.msg(2, 'After block correction, error = %6.2E', norm(P*I - eye(d), 'fro'));
            % one Newton-Raphson iteration
            P = 2*P - P*I*P;
            replab.msg(2, 'After Newton iteration, error = %6.2E', norm(P*I - eye(d), 'fro'));
        end

    end

    methods (Access = protected)

        function sub = computeToSubRep(self)
        % Converts back the generic subrepresentation to a subrepresentation over the original field
            d = self.dimension;
            switch [self.divisionRing '->' self.parent.field]
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
                    [subI, subP] = replab.GenSubRep.correctComplexPair(subI, subP);
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
                    [subI, subP] = replab.GenSubRep.correctQuaternionPair(subI, subP);
                    sub = self.parent.subRep(subI, 'projection', subP, 'divisionAlgebraName', 'quaternion.rep', 'isUnitary', false);
                end
            end
        end

    end


    methods

        function self = GenSubRep(parent, divisionRing, isUnitary, injection, projection)
        % Constructs a generic subrepresentation
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
            sub = self.cached('toSubRep', @() self.computeToSubRep);
        end

        function rho = image(self, g)
            rho = self.projection * self.parent.image(g, 'double/sparse') * self.injection;
        end

        function rho = inverseImage(self, g)
            rho = self.projection * self.parent.inverseImage(g, 'double/sparse') * self.injection;
        end

        function rho = sample(self)
            rho = self.image(self.group.sample);
        end

        function mat = projector(self)
            mat = self.injection * self.projection;
        end

    end

    methods % Implementations

        function l = laws(self)
            l = replab.laws.GenSubRepLaws(self);
        end

    end

    methods (Static)

        function gen = fromQuaternionTypeSubRep(sub)
            assert(sub.overR);
            %assert(strcmp(sub.divisionAlgebraName, 'quaternion.rep'));
            i = replab.H.i;
            j = replab.H.j;
            k = replab.H.k;
            U = [1+i+j+k; 1-i+j-k; 1-i-j+k; 1+i-j-k]/4;
            d = sub.dimension;
            injection = sub.injection('double/sparse') * kron(eye(d/4), U);
            if sub.mapsAreAdjoint
                projection = injection';
                gen = replab.GenSubRep(sub.parent, 'H', true, injection, projection);
            else
                projection = kron(eye(d/4), U)' * sub.projection('double/sparse');
                gen = replab.GenSubRep(sub.parent, 'H', false, injection, projection);
            end
        end

        function gen = fromComplexTypeSubRep(sub)
            assert(sub.overR);
            %assert(strcmp(sub.divisionAlgebraName, 'complex'));
            U = [1; 1i]/sqrt(2);
            d = sub.dimension;
            injection = sub.injection('double/sparse') * kron(eye(d/2), U);
            if sub.mapsAreAdjoint
                projection = injection';
                gen = replab.GenSubRep(sub.parent, 'C', true, injection, projection);
            else
                projection = kron(eye(d/2), U)' * sub.projection('double/sparse');
                gen = replab.GenSubRep(sub.parent, 'C', false, injection, projection);
            end
        end

    end

end
