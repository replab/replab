classdef Declarations < replab.Str
% Describes all the declarations of a class element

    properties
        codeBase % `.CodeBase`: Code base we look in
        classElement % `.ClassElement`: Class element investigated
    end
    
    methods
       
        function self = Declarations(codeBase, classElement)
            self.codeBase = codeBase;
            self.classElement = classElement;
        end
        
        % replab.Str
        
        function str = headerStr(self)
            str = sprintf('All declarations of %s', self.classElement.fullIdentifier);
        end
        
        % Own methods
        
        function els = findAll(self)
        % Returns all declarations of the method/property
            cl = self.classElement.parentClass;
            name = self.classElement.name;
            % concatenate the class under investigation and superclasses
            cls = horzcat({cl}, cl.allSuperclasses);
            % find all classes which have a declaration of our element
            mask = cellfun(@(c) isfield(c.ownElements, name), cls);
            % select those
            cls = cls(mask);
            % find members
            els = cellfun(@(c) c.ownElements.(name), cls, 'uniform', 0);
        end
        
        function els = findDocumentedElements(self)
        % Returns all declarations of the method/property that have documentation 
            els = self.findAll;
            mask = cellfun(@(e) ~e.doc.isempty, els);
            els = els(mask);
        end
        
        function el = findBestDocumented(self)
        % Returns the documented declaration with highest priority
            els = self.findDocumentedElements;
            if isempty(els{1})
                el = [];
            else
                el = els{1};
            end
        end
        
    end

end
