classdef Class < replab.infra.SourceElement

    properties
        superclassIdentifiers
        %ownMethods
        %ownProperties
    end
    
    methods
        
        function self = Class(codeBase, package, classData)
            self = self@replab.infra.SourceElement(codeBase, package, classData.name, 1, ...
                                                   classData.docLines, classData.docLineNumbers);
            sci = {};
            for i = 1:length(classData.superclassIdentifiers)
                id = classData.superclassIdentifiers{i};
                if ~isequal(id, 'handle')
                    sci{1,end+1} = id;
                end
            end
            self.superclassIdentifiers = sci;
        end
        
        function str = name(self)
            str = self.sourceIdentifier;
        end
        
        function p = elementPath(self)
            p = {self.name};
        end

        function asc = allSuperclasses(self)
        % Returns all superclasses of this class, not including itself
            asc = {};
            visitedClassIds = struct;
            queue = self.ownSuperclasses;
            while ~isempty(queue)
                c = queue{1};
                queue = queue(2:end);
                id = c.shmIdentifier;
                if ~isfield(visitedClassIds, id)
                    visitedClassIds.(id) = true;
                    asc{1,end+1} = c;
                    queue = horzcat(queue, c.ownSuperclasses);
                end
            end
        end
        
        function asc = allSubclasses(self)
        % Returns all subclasses of this class, not including itself
            asc = {};
            visitedClassIds = struct;
            queue = self.ownSubclasses;
            while ~isempty(queue)
                c = queue{1};
                queue = queue(2:end);
                id = c.shmIdentifier;
                if ~isfield(visitedClassIds, id)
                    visitedClassIds.(id) = true;
                    asc{1,end+1} = c;
                    queue = horzcat(queue, c.ownSubclasses);
                end
            end
        end
        
        function c = ownSubclasses(self)
            c = self.codeBase.subclasses(self);
        end
        
        function c = ownSuperclasses(self)
            c = cellfun(@(id) self.codeBase.getIdentifier(id), self.superclassIdentifiers, 'uniform', 0);
        end
        
% $$$         
% $$$         function m = member(self, name)
% $$$             m = self.members.(name);
% $$$         end
% $$$         
% $$$         function str = fieldNameEncoding(self)
% $$$             str = strjoin(horzcat(self.packageNameParts, {self.name}), '_');
% $$$         end
% $$$         
% $$$         function elements = findInheritedMember(self, codeBase, name)
% $$$             elements = {};
% $$$             visitedClassIds = {};
% $$$             stack = {self};
% $$$             % breath-first search through the inheritance tree
% $$$             while ~isempty(stack)
% $$$                 % pop item off the stack
% $$$                 h = stack{1};
% $$$                 stack = stack(2:end);
% $$$                 id = h.fieldNameEncoding;
% $$$                 if h.hasMember(name) && ~ismember(id, visitedClassIds)
% $$$                     elements{end+1} = h.member(name);
% $$$                     visitedClassIds{end+1} = id;
% $$$                 end
% $$$                 for i = 1:h.nParents
% $$$                     p = h.parent(codeBase, i);
% $$$                     if ~isempty(p)
% $$$                         stack{end+1} = p;
% $$$                     end
% $$$                 end
% $$$             end
% $$$         end
% $$$         
% $$$         function names = inheritedMemberNames(self, codeBase)
% $$$             names = {};
% $$$             stack = {self};
% $$$             % breath-first search through the inheritance tree
% $$$             while ~isempty(stack)
% $$$                 % pop item off the stack
% $$$                 h = stack{1};
% $$$                 stack = stack(2:end);
% $$$                 % add names of popped item members
% $$$                 names = union(names, h.memberNames);
% $$$                 % add parents to the stack
% $$$                 for i = 1:h.nParents
% $$$                     p = h.parent(codeBase, i);
% $$$                     if ~isempty(p)
% $$$                         stack{end+1} = p;
% $$$                     end
% $$$                 end
% $$$             end
% $$$         end
% $$$         

    end
    
end
