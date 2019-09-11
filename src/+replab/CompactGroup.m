classdef CompactGroup < replab.Group
% A group equipped with a Haar measure

    methods
        
        %% Abstract
        
        function g = sampleUniformly(self)
        % Samples from the Haar measure
        %
        % Returns:
        %   element: Group element sampled from the Haar measure
            error('Abstract');
        end
        
        %% Representations
        
        function rep = trivialRep(self, field, dimension)
        % Returns the trivial representation of this group on a finite dimensional vector space
        %
        % For convenience, either the representation can act on a real or complex vector space,
        % and multiple copies of the 1-dimensional trivial representation can be included, when
        % dimension > 1.
        % 
        % Args:
        %   field ({'R', 'C'}): Whether the representation is real (R) or complex (C)
        %   dimension (integer): Representation dimension
        %
        % Returns:
        %   replab.Rep: An instance of the trivial representation
            rep = replab.rep.TrivialRep(self, field, dimension);
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
