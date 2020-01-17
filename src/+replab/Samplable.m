classdef Samplable < replab.Str
% Describes a set of elements that can be sample
%
% Random elements are produced by the method `.sample`.

    methods % ABSTRACT

        function t = sample(self)
        % Samples an element from this domain
        %
        % This method does not make any guarantees about genericity, and is primarily used for law checks.
        %
        % Returns:
        %   element: Random domain element
            error('Abstract');
        end

    end

    methods (Static)

        function s = lambda(header, sampleFun)
        % Constructs a domain from function handles
        %
        % Args:
        %   header (char): Header display string
        %   eqvFun (function_handle): Handle implementing the `eqv` method
        %   sampleFun (function_handle): Handle implementing the `sample` method
        %
        % Returns:
        %   `+replab.Samplable`: The constructed domain
            s = replab.lambda.Samplable(header, sampleFun);
        end

    end

end
