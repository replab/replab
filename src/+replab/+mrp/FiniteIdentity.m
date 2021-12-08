classdef FiniteIdentity < replab.FiniteIsomorphism & replab.mrp.Identity
% Identity isomorphism

    methods

        function self = FiniteIdentity(group)
            self@replab.mrp.Identity(group);
        end

    end

    methods % Implementations

        % Obj

        function l = laws(self)
            l = laws@replab.FiniteIsomorphism(self);
        end
        
        % Morphism

        function y = imageElement(self, x)
            y = x;
        end

        % Isomorphism

        function x = preimageElement(self, y)
            x = y;
        end

        % FiniteMorphism

        function T = imageGroup(self, S)
            T = S;
        end

        function S = preimageGroup(self, T)
            S = T;
        end

        % FiniteIsomorphism

        function l = preservesTypeOrder(self)
            l = true;
        end

        % Str

        function [names, values] = additionalFields(self)
            [names, values] = additionalFields@replab.FiniteIsomorphism(self);
        end

    end

    methods (Access = protected) % Implementations

        function inv = computeInverse(self)
            inv = self;
        end

    end

    methods (Static)

        function m = identity(group)
            m = identity@replab.FiniteIsomorphism(group);
        end
        
        function m = lambda(source, target, preimageElementFun, imageElementFun, torusMap)
            m = lambda@replab.FiniteIsomorphism(source, target, preimageElementFun, imageElementFun, torusMap);
        end

    end

end
