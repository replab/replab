classdef Equivariant < replab.Domain
% Describes the vector space of group-equivariant maps
%
% Let repR and repC be two representations of the same group G.
% This describes the set of matrices X such that
%
% repR.image(g) * X = X * repC.image(g)
%   
% See Proposition 4 of 
% J.-P. Serre, Linear Representations of Finite Groups (Springer, 1977).
    
    properties (SetAccess = protected)
        field % vector space field, 'R' (real) or 'C' (complex)
        nR    % row size
        nC    % column size
        group % group
        repR  % representation of row space
        repC  % representation of column space
    end
    
    properties (Access = protected)
        parent_ % parent domain: real or complex matrices
    end
    
    methods
       
        function self = Equivariant(repR, repC)
            
        end
    end
         
end