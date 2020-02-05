classdef TrivialInfo < replab.irreducible.Info

    methods

        function self = TrivialInfo(field)
            switch field
              case 'R'
                da = 'R';
              case 'C'
                da = [];
              otherwise
                error('Unknown field %s', field);
            end
            self@replab.irreducible.Info(da, []);
        end

    end

end
