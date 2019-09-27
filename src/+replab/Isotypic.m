classdef Isotypic < replab.Str
% Describes an isotypic component in the decomposition of a representation
%
% An isotypic component regroups equivalent irreducible representations, expressed in the same basis.
% Note that if the multiplicity is not one, there is a degeneracy in the basis of the copies, and
% the particular basis choosen is not deterministic.
%
% However the subspace spanned by an isotypic component as a whole is unique.
    
    properties
        parent % replab.Rep: Representation of which this is an isotypic component
        copies % row cell array of replab.SubRep: Equivalent irreducible subrepresentations in
               %                                  this isotypic component
        multiplicity % integer: number of copies in this isotypic component
        copyDimension % integer: dimension of each irreducible representation in this component
    end
    
    methods
        
        function self = Isotypic(parent, copies)
            assert(isa(parent, 'replab.Rep'));
            assert(length(copies) >= 1, 'Isotypic component cannot be empty');
            for i = 1:length(copies)
                ci = copies{i};
                assert(isa(ci, 'replab.SubRep'));
                assert(ci.isKnownCanonicalIrreducible);
            end
            self.parent = parent;
            self.copies = copies;
            self.multiplicity = length(copies);
            self.copyDimension = copies{1}.dimension;
        end
        
        function d = dimension(self)
        % Returns the dimension of this isotypic component
            d = self.copyDimension * self.multiplicity;
        end
        
        function r = rep(self)
        % Returns the subrepresentation corresponding to this isotypic component
            U = zeros(0, self.parent.dimension);
            for i = 1:self.nCopies
                U = [U;self.copy(i).U];
            end
            r = self.parent.subRepUnitary(U);
            % TODO: preserve rational bases
        end
        
        function n = nCopies(self)
        % Returns the number of copies = the multiplicity
            n = self.multiplicity;
        end
        
        function c = copy(self, i)
        % Returns the i-th copy of the irreducible representation
            c = self.copies{i};
        end
        
        function names = hiddenFields(self)
            names = hiddenFields@replab.Str(self);
            names{1, end+1} = 'copies';
        end
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Str(self);
            for i = 1:self.nCopies
                names{1, end+1} = sprintf('copy(%d)', i);
                values{1, end+1} = self.copy(i);
            end
        end
        
        function s = headerStr(self)
            if isequal(self.copy(1).field, 'C')
                rt = 'C';
            else
                rt = self.copy(1).irrepInfo.divisionAlgebra;
                assert(~isempty(rt));
            end
            if self.multiplicity > 1
                s = sprintf('Isotypic component I(%d)x%s(%d)', self.multiplicity, rt, self.copyDimension);
            else
                s = sprintf('Isotypic component %s(%d)', rt, self.copyDimension);
            end
        end
        
    end
    
end
