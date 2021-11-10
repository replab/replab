classdef NiceIsomorphism < replab.FiniteIsomorphism
% An isomorphism from a subgroup of a finite group type to a finite group (Abstract)
%
% In addition, the isomorphism should preserve the type order.

    methods (Static)

        function targetArgs = translateKeywordArgs(sourceArgs)
        % Translates the key/value pairs given for a source group for the target (nice) group
        %
        % Args:
        %   sourceArgs (cell(1,\*)): Key/value pairs
        %
        % Returns:
        %   cell(1,\*): Key/value pairs
            targetArgs = cell(1, 0);
            i = 1;
            while i < length(sourceArgs)
                switch sourceArgs{i}
                  case 'generatorNames'
                    targetArgs{1,end+1} = 'generatorNames';
                    targetArgs{1,end+1} = sourceArgs{i+1};
                  case 'order'
                    targetArgs{1,end+1} = 'order';
                    targetArgs{1,end+1} = sourceArgs{i+1};
                  case 'relators'
                    targetArgs{1,end+1} = 'relators';
                    targetArgs{1,end+1} = sourceArgs{i+1};
                  case 'abelianInvariants'
                    targetArgs{1,end+1} = 'abelianInvariants';
                    targetArgs{1,end+1} = sourceArgs{i+1};
                    % TODO: conjugacyClasses
                end
                i = i + 2;
            end
        end

    end

    methods (Access = protected)

        function finishConstruction(self, sourceFun, sourceGenerators, targetType)
        % Finishes the construction
        %
        % When ``sourceFun`` is called, note that this object does not have the `.source` property set.
        %
        % Args:
        %   sourceFun (function_handle): Function handle that takes this object as argument and returns the source group
        %   targetType (`+replab.FiniteGroupType`): Type of target finite group
            targetGenerators = cellfun(@(g) self.imageElement(g), sourceGenerators, 'uniform', 0);
            self.target = targetType.groupWithGenerators(targetGenerators);
            self.torusMap = [];
            self.source = sourceFun(self);
        end

    end

    methods % Implementations

        % Obj

        function l = laws(self)
            l = replab.laws.OrderPreservingFiniteIsomorphismLaws(self);
        end

    end


    methods

        function t = sourceType(self)
        % Returns the type of the finite group which is the source of this morphism
            t = self.source.type;
        end

        function l = sourceContains(self, s)
        % Returns whether the source of this morphism contains the given type element
        %
        % Args:
        %   s (element of `.sourceType`): Element to check
        %
        % Returns:
        %   logical: True if source contains the given element
            error('Abstract');
        end

    end

end
