classdef Domain < replab.Samplable
% Describes a set of elements with a common structure
%
% At the base of the hierarchy, `~replab.Domain` describes a set of elements
% that can be tested for equality (`~replab.Domain.eqv`) and from which random samples
% can be taken (`~replab.Domain.sample`). Such sets are potentially infinite.
%
% As `~replab.Domain` is an abstract base class, it contains abstract methods.
%
% To quickly create an instance of `~replab.Domain`, the method `~replab.Domain.lambda` can be used,
% passing the method implementations as function handles.

% Those elements can be compared (`~replab.Domain.eqv`), and random elements can be produced (`~replab.Domain.sample`).

    methods % ABSTRACT

        function b = eqv(self, t, u)
        % Tests domain elements for equality/equivalence
        %
        % Args:
        %   t (domain element): First element to test
        %   u (domain element): Second element to test
        %
        % Returns:
        %   logical: True when ``t`` and ``u`` are equivalent, and false otherwise
            error('Abstract');
        end

    end

    methods

        function h = hash(self, x)
        % Returns a hash code value for the given domain element
        %
        % The value returned is of double type, and two cases are possible:
        %
        % - a valid hash code can be computed for the object, and the method returns
        %   a nonnegative integer between 0 and 2^32-1
        % - no valid hash code can be computed, in which case the method returns NaN
        %
        %
        % Args:
        %   x (domain element): Element to hash
        %
        % Returns:
        %   double: Hash value as a double encoding an unsigned 32 bits integer, or NaN
            h = NaN;
        end

    end

    methods (Static)

        function domain = lambda(header, eqvFun, sampleFun)
        % Constructs a domain from function handles
        %
        % Args:
        %   header (char): Header display string
        %   eqvFun (function_handle): Handle implementing the `eqv` method
        %   sampleFun (function_handle): Handle implementing the `sample` method
        %
        % Returns:
        %   `+replab.Domain`: The constructed domain
            domain = replab.lambda.Domain(header, eqvFun, sampleFun);
        end

        function h = hashVector(v)
        % Implements the Jenkins one-at-the-time hash function
        %
        % Args:
        %   v (double(1,\*)): Row vector of values to hash, the values must be between 0 and 2^32-1, or be NaN
        %
        % Returns:
        %   double: NaN is any element of v is NaN, otherwise an integer between 0 and 2^32-1 encoded as a double
            if any(isnan(v))
                h = NaN;
            else
                mask = uint64(2^32-1);
                h = uint64(0);
                for i = 1:length(v)
                    h = h + uint64(v(i));
                    h = h + bitshift(h, 10);
                    h = bitand(h, mask);
                    h = bitxor(h, bitshift(h, -6));
                end
                h = h + bitshift(h, 3);
                h = bitxor(h, bitshift(h, -11));
                h = h + bitshift(h, 15);
                h = double(bitand(h, mask));
            end
        end

    end

    methods % Test helpers

        function assertNotEqv(self, x, y, context)
        % Compares two elements for inequality
        %
        % Asserts that ``x`` and ``y`` are equivalent
        %
        % Args:
        %   x (domain element): first element
        %   y (domain element): second element to compare
        %   context (charstring): context
            if self.eqv(x, y)
                errorDesc = 'The values %s and %s are equivalent, but they should not be';
                errorId = 'assertNotEqual:equal';
            else
                return
            end

            if nargin < 4
                context = '';
            end

            names = evalin('caller', 'who');
            nV = length(names);
            values = cell(1, nV);
            for i = 1:nV
                values{i} = evalin('caller', names{i});
            end

            message = replab.laws.message(errorDesc, context, {x y}, names, values);
            if moxunit_util_platform_is_octave()
                error(errorId, '%s', message);
            else
                throwAsCaller(MException(errorId, '%s', message));
            end

        end

        function assertEqv(self, x, y, context)
        % Compares two elements for equality
        %
        % Asserts that ``x`` and ``y`` are not equivalent
        %
        % Args:
        %   x (domain element): first element
        %   y (domain element): second element to compare
        %   context (charstring, optional): context
            if ~self.eqv(x, y)
                errorDesc = 'The values %s and %s are not equivalent';
                errorId = 'assertEqual:nonEqual';
            else
                return
            end

            if nargin < 4
                context = '';
            end

            names = evalin('caller', 'who');
            nV = length(names);
            values = cell(1, nV);
            for i = 1:nV
                values{i} = evalin('caller', names{i});
            end

            message = replab.laws.message(errorDesc, context, {x y}, names, values);
            if moxunit_util_platform_is_octave()
                error(errorId, '%s', message);
            else
                throwAsCaller(MException(errorId, '%s', message));
            end

        end

    end

end
