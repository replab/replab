classdef Samplable < replab.Obj
% Describes a set of elements that can be sample
%
% Random elements are produced by the method `.sample`.

    methods % ABSTRACT

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
