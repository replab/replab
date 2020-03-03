classdef Group < testpkg.Monoid
% Describes a group
    
    methods % Abstract methods
        
        function xInv = inverse(self, x)
        % Computes the inverse of an element
        %
        % Given ``x``, returns ``xInv`` such that
        %
        % ``x xInv = identity``
        % 
        % Args:
        %   x (element): Group element to compute the inverse of
        %
        % Returns:
        %   element: Inverse of ``x``
            error('Abstract');
        end
        
    end

end
