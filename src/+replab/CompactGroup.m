classdef CompactGroup < replab.Group
% A group equipped with a Haar measure

    methods
        
        function g = sampleUniformly(self)
        % Samples from the Haar measure
        %
        % Returns:
        %   element: Group element sampled from the Haar measure
            error('Abstract');
        end

    end

    methods (Static)
        
        function group = lambda(header, eqvFun, sampleFun, ...
                                composeFun, identity, inverseFun, ...
                                sampleUniformlyFun)
        % Constructs a compact group from function handles
        %
        % Args:
        %   header (char): Header display string
        %   eqvFun (function_handle): Handle implementing the `eqv` method
        %   sampleFun (function_handle): Handle implementing the `sample` method
        %   composeFun (function_handle): Handle implementing the `compose` method
        %   identity (element): Identity element of this monoid
        %   inverseFun (function_handle): Handle implementing the `inverse` method
        %   sampleUniformlyFun (function_handle): Handle implementing the `sampleUniformly` method
        %
        % Returns:
        %   replab.CompactGroup: The constructed compact group
            
            group = replab.lambda.CompactGroup(header, eqvFun, sampleFun, ...
                                               composeFun, identity, inverseFun, ...
                                               sampleUniformlyFun);
        end
        
    end
    
end
