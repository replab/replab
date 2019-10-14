classdef Group < replab.Monoid
% Describes a group
    
    methods % Abstract methods
        
        function xInv = inverse(self, x)
        % Computes the inverse of an element
        %
        % Given `x`, returns `xInv` such that
        %
        % `x xInv = identity`
        % 
        % Args:
        %   x (element): Group element to compute the inverse of
        %
        % Returns:
        %   element: Inverse of `x`
            error('Abstract');
        end
        
    end

    methods % Methods with default implementations
        
        function x = leftConjugate(self, by, on)
        % Returns the left conjugate of a group element
        %
        % Convenience method that can be overriden for speed optimizations
        %
        % Args:
        %   by (element): Element conjugating
        %   on (element): Element conjugated
        %  
        % Returns:
        %   element: left conjuagate, i.e. `by * on * by^-1` in multiplicative notation
            x = self.composeWithInverse(self.compose(by, on), by);
        end
        
        function z = composeWithInverse(self, x, y)
        % Returns the composition of an element with the inverse of another element
        %
        % Convenience method that can be overriden for speed optimizations
        %
        % Args:
        %   x (element): First element
        %   y (element): Second element
        %
        % Returns:
        %   element: the result of ``x * y^-1`` in multiplicative notation
            z = self.compose(x, self.inverse(y));
        end
        
    end

    methods (Static)
        
        function group = lambda(header, eqvFun, sampleFun, ...
                                composeFun, identity, inverseFun)
        % Constructs a group from function handles
        %
        % Args:
        %   header (char): Header display string
        %   eqvFun (function_handle): Handle implementing the `eqv` method
        %   sampleFun (function_handle): Handle implementing the `sample` method
        %   composeFun (function_handle): Handle implementing the `compose` method
        %   identity (element): Identity element of this monoid
        %   inverseFun (function_handle): Handle implementing the `inverse` method
        %
        % Returns:
        %   :class:`+replab.Group`: The constructed group
            
            group = replab.lambda.Group(header, eqvFun, sampleFun, ...
                                        composeFun, identity, inverseFun);
        end
        
    end
    
end
