classdef Samplable < replab.Obj
% Describes a set of elements that can be sample
%
% Random elements are produced by the method `.sample`.

    methods % ABSTRACT

        function t = sample(self)
        % Samples an element from this set
        %
        % In general, this method does not make guarantees about genericity.
        %
        % For `.CompactGroup` however, this method must sample uniformly from the Haar measure.
        %
        % Returns:
        %   set element: Random set element
            error('Abstract');
        end

    end

    methods (Static)

        function s = lambda(header, sampleFun)
        % Constructs a samplable set from function handles
        %
        % Args:
        %   header (char): Header display string
        %   sampleFun (function_handle): Handle implementing the `.sample` method
        %
        % Returns:
        %   `+replab.Samplable`: The constructed samplable set
            s = replab.lambda.Samplable(header, sampleFun);
        end

    end

end
