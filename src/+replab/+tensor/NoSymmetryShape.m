classdef NoSymmetryShape < replab.tensor.Shape
% Describes the shape of a column-major ordered tensor without symmetries
    
    properties (Access = protected)
        cumProd % uint32: equal to the first `.rank` elements of cumprod(self.dimensions)
    end
    
    methods
        
        function self = NoSymmetryShape(dimensions, group, isOrderColumnMajor)
            assert(group.isTrivial);
            assert(isOrderColumnMajor);
            self = self@replab.tensor.Shape(dimensions, group, isOrderColumnMajor);
            cp = cumprod(uint32(dimensions));
            self.cumProd = [uint32(1) cp(1:end-1)];
       end
       
       function d = nComponents(self)
           d = prod(self.dimensions);
       end

       function sub = indToSub(self, ind)
           ind = ind(:);
            nRows = length(ind);
            switch self.rank
              case 0
                sub = zeros(nRows, 0, class(ind));
              case 1
                sub = ind;
              otherwise
                indm = uint32(ind - 1); % convert to 0-based indices
                dims = uint32(self.dimensions);
                sub = zeros(nRows, self.rank, 'uint32');
                for i = 1:self.rank
                    sub(:, i) = mod(indm, dims(i));
                    indm = indm - sub(:, i);
                    indm = indm ./ dims(i);
                end
                sub = cast(sub + 1, class(ind)); % convert back to 1-based indices
            end
        end
        
        function ind = subToInd(self, sub)
            nRows = size(sub, 1);
            switch self.rank
              case 0
                ind = ones(nRows, 1, class(sub));
              case 1
                ind = sub;
              otherwise
                assert(size(sub, 2) == self.rank);
                if isa(sub, 'uint32')
                    ind = zeros(nRows, 1, 'uint32');
                    ind = sub(:, 1);
                    for i = 2:n
                        ind = ind + self.cumProd(i) * (sub(:, i) - 1);
                    end
                else
                    ind = cast((double(self.cumProd) * (sub' - 1) + 1)', class(sub));
                end
            end
        end
        
        function orbit = subOrbit(self, sub)
            orbit = sub; % only one element in orbit
        end
        
   end
   
end
