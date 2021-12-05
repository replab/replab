classdef DirectProductGroup_finite < replab.DirectProductGroup & replab.gen.FiniteGroup
% External direct product of finite groups
%
% In particular, the permutation image of an element of a direct product group
% is simply the concatenation of the permutation images of the factors (which
% are nice finite groups themselves).
%
% We overload a bunch of methods to make sure we use the more efficient variants, that do not require the BSGS chain construction.

    methods

        function self = DirectProductGroup_finite(factors, type, generators, nice, niceIsomorphism)
            self@replab.gen.FiniteGroup(type, generators, 'nice', nice, 'niceIsomorphism', niceIsomorphism);
            self.factors = factors;
        end

    end

    methods % Implementations

        % Str

        function s = headerStr(self)
            s = headerStr@replab.DirectProductGroup(self);
        end

        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.DirectProductGroup(self), ...
                hiddenFields@replab.gen.FiniteGroup(self) ...
                );
        end

        function [names values] = additionalFields(self)
            [names1 values1] = additionalFields@replab.DirectProductGroup(self);
            [names2 values2] = additionalFields@replab.gen.FiniteGroup(self);
            names = replab.str.horzcatForce(names1, names2);
            values = replab.str.horzcatForce(values1, values2);
        end

        % Domain

        function b = eqv(self, x, y)
            b = eqv@replab.DirectProductGroup(self, x, y);
        end

        function g = sample(self)
            g = sample@replab.DirectProductGroup(self);
        end

        % Monoid

        function z = compose(self, x, y)
            z = compose@replab.DirectProductGroup(self, x, y);
        end

        % Group

        function xInv = inverse(self, x)
            xInv = inverse@replab.DirectProductGroup(self, x);
        end

        % FiniteSet

        function s = setProduct(self)
            T = {};
            % The decomposition of a direct product into sets
            % is simply the concatenation of the sequence of sets
            % corresponding to each factor
            for i = 1:self.nFactors
                D = self.factor(i).setProduct.sets;
                Ti = cell(1, length(D));
                for j = 1:length(D)
                    Dj = D{j};
                    Tij = cell(1, length(Dj));
                    for k = 1:length(Dj)
                        Djk = Dj{k};
                        Tijk = self.identity;
                        Tijk{i} = Djk;
                        Tij{k} = Tijk;
                    end
                    Ti{j} = Tij;
                end
                T = horzcat(T, Ti);
            end
            s = replab.SetProduct(self, T, true);
        end

    end


% $$$     methods (Access = protected) % Implementations
% $$$
% $$$         function g = atFun(self, ind)
% $$$        % See comments in self.elementsSequence
% $$$             g = self.identity;
% $$$             ind = ind - 1;
% $$$             for i = self.nFactors:-1:1
% $$$                 f = self.factor(i);
% $$$                 this = mod(ind, f.order);
% $$$                 ind = (ind - this)/f.order;
% $$$                 g{i} = f.elementsSequence.at(this + 1);
% $$$             end
% $$$         end
% $$$
% $$$         function ind = findFun(self, g)
% $$$         % See comments in self.elementsSequence
% $$$             ind = vpi(0);
% $$$             for i = 1:self.nFactors
% $$$                 f = self.factor(i);
% $$$                 ind = ind * f.order;
% $$$                 ind = ind + f.elementsSequence.find(g{i}) - 1;
% $$$             end
% $$$             ind = ind + 1;
% $$$         end
% $$$
% $$$         function c = computeCharacterTable(self, field)
% $$$             c = replab.ct.directProduct(self, field);
% $$$         end
% $$$
% $$$         function c = computeConjugacyClasses(self)
% $$$             classes = cellfun(@(f) f.conjugacyClasses.classRepresentatives, self.factors, 'uniform', 0);
% $$$             reps = replab.util.cartesian(classes);
% $$$             c = replab.ConjugacyClasses(self, cellfun(@(r) self.conjugacyClass(r), reps, 'uniform', 0));
% $$$         end
% $$$
% $$$         function e = computeElementsSequence(self)
% $$$             e = replab.Sequence.lambda(self.order, ...
% $$$                                        @(ind) self.atFun(ind), ...
% $$$                                        @(g) self.findFun(g));
% $$$             % The elements of a direct product of finite groups can be
% $$$             % enumerated by considering the direct product as a cartesian
% $$$             % product of sets, and decomposing the index a la ind2sub/sub2ind
% $$$             % which is the role of the `atFun` and `findFun` functions
% $$$
% $$$             % TODO: verify ordering
% $$$         end
% $$$
% $$$
% $$$     end

% $$$     methods % Implementations
% $$$
% $$$
% $$$     end

end
