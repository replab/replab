classdef CodeBase < replab.Str
    
    properties
        packages % struct-based hash map
    end
    
    methods
        
        function p = lookupPackage(self, nameParts)
        % Looks for a package from its name parts
        %
        % Args:
        %   nameParts (cell row vector of charstring): Parts of the package name
        %
        % Returns:
        %   :class:`+replab.+infra.Package`: The corresponding package, or ``[]`` if not found
            fname = strjoin(nameParts, '_');
            if isfield(self.packages, fname)
                p = getfield(self.packages, fname);
            else
                p = [];
            end
        end
        
        function [package packageNameParts restNameParts] = lookupPackageGreedy(self, nameParts)
        % Greedily looks up for a package from its name parts, disregarding the suffix that does not match 
        %
        % Returns `package` which matches `packageNameParts`, and `restNameParts` such that
        % ``nameParts = horzcat(packageNameParts, restNameParts)``.
        %
        % Args:
        %   nameParts (cell row vector of charstring): Parts of a fully qualified object name
        %
        % Returns
        % -------
        %   package: 
        %     :class:`+replab.+infra.Package`: The package looked up for
        %   packageNameParts:
        %     cell row vector of charstring: The maximal prefix of `nameParts` that matches a package name
        %   restNameParts:
        %     cell row vector of charstring: The tail of `nameParts` that could not be matched
            n = length(nameParts);
            for i = n:-1:0
                packageNameParts = nameParts(1:i);
                restNameParts = nameParts(i+1:n);
                package = self.lookupPackage(packageNameParts);
                if ~isempty(package)
                    return
                end
            end
            error('Should not happen: empty name parts match the root package');
        end

        function obj = lookupName(self, name)
        % Finds an object from its fully qualified name
        %
        % Args:
        %   name (charstring): Name of the object to lookup
        %
        % Returns:
        %   The looked up object
        %
        % Raises:
        %   An error if the object cannot be found or the name is malformed
            nameParts = strsplit(name, '.');
            [package packageNameParts restNameParts] = self.lookupPackageGreedy(nameParts);
            obj = package.lookupMemberName(restNameParts);
        end
        
    end
    
end
