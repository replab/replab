classdef Shape
% Describes the dimensions of a tensor and its possible symmetries
%
% It provides efficient implementations of ind2sub/sub2ind vectorizing operations that
% avoid the use of cell arrays, and take in account symmetry to reduce storage.
%
% The tensor elements can either be stored in column-major or row-major order. See the
% `Wikipedia article <https://en.wikipedia.org/wiki/Row-_and_column-major_order>`_ for
% details.
%
% In particular, MATLAB/Octave uses column major ordering, as does Fortran; while row major
% ordering is used in C/C++; row major ordering also respect the ordering of the Kronecker
% product.
%
% Use `+replab.+tensor.Shape.make` to construct a `Shape` object.
    
    properties (SetAccess = immutable)
        dimensions % integer row vector: Dimensions of the tensor
        group % `+replab.PermutationGroup`: Index symmetry group
        isOrderColumnMajor % logical: Whether the linear ordering is column major
    end
    
    properties (Access = protected)
        unsymmetrized_ = [];
    end
    
    methods (Access = protected)
        
        function self = Shape(dimensions, group, isOrderColumnMajor)
            assert(isvector(dimensions));
            assert(isa(dimensions, 'double'));
            assert(isa(group, 'replab.PermutationGroup'));
            assert(group.domainSize == length(dimensions));
            self.dimensions = dimensions;
            self.group = group;
            self.isOrderColumnMajor = isOrderColumnMajor;
        end

    end
    
    methods
        
        function n = rank(self)
        % Returns the number of indices in the tensor
        %
        % ``n = shape.rank``
        %
        % Returns:
        %   integer: The tensor rank
            n = length(self.dimensions);
        end
        
        function n = nComponents(self)
        % Returns the number of independent components in this tensor
        %
        % ``n = shape.nComponents``
        %
        % Returns:
        %   integer: The number of independent components
            error('Abstract');
        end
        
        function shape = unsymmetrized(self)
        % Returns this shape without any of its symmetries
        %
        % ``unsym = shape.unsymmetrized``
        %
        % Returns:
        %   `+replab.+tensor.Shape`: The shape with the same dimensions and ordering type, but with ``group.order = 1``
            if isempty(self.unsymmetrized_)
                group = self.group.trivialSubgroup;
                self.unsymmetrized_ = replab.tensor.Shape.make(self.dimensions, group, self.isOrderColumnMajor);
            end
            shape = self.unsymmetrized_;
        end
        
        function can = canonicalSub(self, sub)
        % Returns the canonical form of subindices
        %
        % See `.indToSub` regarding argument and return types.
        %
        % Args:
        %   sub (integer matrix): A ``m x self.rank`` integer matrix of subindices
        %
        % Returns:
        %   integer matrix: A ``m x self.rank`` integer matrix of canonical subindices
            can = self.indToSub(self.subToInd(sub));
        end
        
        function sub = indToSub(self, ind)
        % Converts from unique component indices to subindices
        %
        % The internal implementation can use different types; the input argument can be
        % either of type ``double`` or of type ``uint32`` (the latter being usually more
        % efficient).
        % The return argument type matches the input argument type, and the returned
        % subindices are canonical.
        %
        % Args:
        %   ind (integer column vector): A vector of ``m`` indices with values in ``1...self.nComponents``
        %
        % Returns:
        %   integer matrix: A ``m x self.rank`` integer matrix of subindices
            error('Abstract');
        end
        
        function ind = subToInd(self, sub)
        % Converts from subindices to unique component indices
        %
        % See `.indToSub` regarding argument and return types.
        % 
        % Args:
        %   sub (intger matrix): A ``m x self.rank`` integer
        %
        % Returns:
        %   integer column vector: A vector of ``m`` indices with values in ``1...self.nComponents``
            error('Abstract');
        end
        
        function orbit = subOrbit(self, sub)
        % Returns all subindices equivalent to the given subindices
        %
        % Args:
        %   sub (integer row vector): A vector of ``self.rank`` indices
        %
        % Returns:
        %   integer matrix: A matrix whose rows are the orbit elements
            error('Abstract');
        end

    end
    
    methods (Static)
       
        function shape = make(dimensions, group, isOrderColumnMajor)
            if ~isOrderColumnMajor
                shape = replab.tensor.RowMajorShape(dimensions, group, isOrderColumnMajor);
                return
            end
            partition = group.orbits;
            % Now, we assume that isOrderColumnMajor == true
            if group.order == replab.Permutations(length(dimensions)).order
                shape = replab.tensor.FullSymmetryShape(dimensions, group, isOrderColumnMajor);
            elseif group.isTrivial
                shape = replab.tensor.NoSymmetryShape(dimensions, group, isOrderColumnMajor);
            end
        end
        
    end
    
end
