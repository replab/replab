classdef SourceElement < replab.infra.Element

    properties
        package % `replab.infra.Package`: Package this element is part of
        sourceIdentifier % charstring: Name of the ``.m`` file without extension
        startLineNumber % integer: Line number at which this object starts in the source file
        doc % `.Doc`: Documentation comment
    end
    
    methods
        
        function self = SourceElement(codeBase, package, sourceIdentifier, startLineNumber, docLines, docLineNumbers)
            self = self@replab.infra.Element(codeBase);
            self.package = package;
            self.sourceIdentifier = sourceIdentifier;
            self.startLineNumber = startLineNumber;
            self.doc = replab.infra.Doc(self, docLines, docLineNumbers);
        end
        
        function p = elementPath(self)
            error('Abstract');
        end
        
        function [packagePath elementPath] = splitPath(self)
            packagePath = self.package.packagePath;
            elementPath = self.elementPath;
        end
        
        function str = sourceLinkOpen(self)
            str = replab.infra.linkOpen('%s:%d', '%s:%d', self.filename, self.startLineNumber);
        end
        
        function fn = absoluteFilename(self)
        % Returns the full path to the file containing the source code of this object
            parts = self.relativeFilenameParts;
            fn = fullfile(self.codeBase.rootFolder, parts{:});
        end
        
        function parts = relativeFilenameParts(self)
        % Returns the relative path to the source code of this object
        %
        % For example, this could return ``{'+replab' 'Group.m'}``
            parts = horzcat(cellfun(@(x) ['+' x], self.packagePathParts, 'uniform', 0), {[self.sourceIdentifier '.m']});
        end

    end

end
