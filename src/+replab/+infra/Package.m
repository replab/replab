classdef Package < replab.Str
    
    properties
        nameParts % row cell vector of string: parts of the package
        members % struct-based hash map
    end

    methods
       
        function self = Package(nameParts, members)
            self.nameParts = nameParts;
            self.members = members;
        end

        function n = nMembers(self)
            n = length(self.members);
        end
        
        function m = member(self, i)
            m = self.members{i};
        end
        
        function names = hiddenFields(self)
            names = hiddenFields@replab.Str(self);
            names{1, end+1} = 'members';
        end
        
        function [names values] = additionalFields(self)
            [names values] = additionalFields@replab.Str(self);
            for i = 1:self.nMembers
                names{1, end+1} = sprintf('member(%d)', i);
                values{1, end+1} = self.member(i);
            end
        end

        function s = headerStr(self)
            s = ['Package ' strjoin(self.nameParts, '.')];
        end
        
        function member = lookupMemberName(self, nameParts)
        % Looks up a package member from its name parts
        %
        % Note that subpackages cannot be looked up by this method.
        %
        % Args:
        %   nameParts (cell row vector of charstring): Parts of the member name
        %
        % Returns:
        %   The member corresponding to the name parts
        %   
        % Raises:
        %   An error if the object cannot be found or the name is malformed
            if isempty(nameParts)
                member = self;
            else
                elementName = nameParts{1};
                if isfield(self.elements, elementName)
                    element = getfield(self.members, elementName);
                    if isa(element, 'replab.infra.Function')
                        if length(nameParts) > 1
                            error('A function does not have members');
                        end
                        member = element;
                        return
                    elseif isa(element, 'replab.infra.Class')
                        switch length(nameParts)
                          case 1
                            member = element;
                          case 2
                            member = element.lookupMemberName(nameParts{2});
                          otherwise
                            error('Members of a class do not have submembers');
                        end
                    end
                end
            end
        end
        
    end
    
end
