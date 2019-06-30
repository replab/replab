function X = prod(A, dim)
% PROD   Product of array elements.
% (Quaternion overloading of standard Matlab function.)
%
% The product is computed with the low index elements on the left and the
% high index elements on the right.

% Copyright (c) 2008, 2015 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 2), nargoutchk(0, 1) 

if nargin == 1

    % There is no dim parameter, so if the input is a vector, we can just
    % multiply all the elements together, if it is a matrix or array, we
    % need to work on the columns,....
    
    if isvector(A)
        X = product(A);
    elseif ismatrix(A)
        % The result will be a row vector, one entry per column of A, so we
        % make a row vector from the first row of A and replace its entries
        % with the products down the columns. We need to use subscripted
        % references, but these do not work inside a class method (which
        % this is). See the file 'implementation notes.txt', item 8.
        
        X = subsref(A, substruct('()', {1, ':'})); % X = A(1, :); 
        
        for k = 1:length(X)
            % X(k) = product(A(:, k));
            X = subsasgn(X, substruct('()', {k}), ...
                product(subsref(A, substruct('()', {':', k}))));
        end
        
    else % A has more than 2 dimensions.
        
        % Find the first non-singleton dimension and make a recursive call
        % with an explicit dimension parameter.
        
        d = size(A);
        for dim = 1:length(d)
            if d(dim) ~= 1
                X = prod(A, dim);
                return;
            end
        end
        
        % This should not happen, since if none of the dimensions is > 1,
        % the array should not have ndims > 2. But just in case ......
        
        error('There are no non-singleton dimensions to the array!')
    end
    
else
    if ~isnumeric(dim)
        error('Dimension argument must be numeric');
    end
    
    if ~isscalar(dim) || ~ismember(dim, 1:ndims(A))
        error(['Dimension argument must be a positive'...
               ' integer scalar within indexing range.']);
    end
    
    if ismatrix(A)
        
        % A is a matrix and it can be handled by the code above. There are
        % two cases since the dim parameter can be only 1 or 2.
        
        if dim == 1
            % We are required to compute along the columns, which is the
            % default, so a plain recursive call will do the trick.
            
            X = prod(A);
        else
        
            % We are required to compute the product along the rows and
            % return a column vector. We can handle this using transpose, a
            % recursive call to the code above, and transposition of the
            % result.

            X = prod(A.').';
        end
        return
    end
    
    % TODO: The more general case of multiple dimensions is more difficult.
    % The problem is how to access one dimension of the array when the
    % number of dimensions is variable. Possible approaches might be to use
    % permute, reshape, shiftdims etc. But how do you access one dimension?
    % You cannot write A(1,2,:,3,4) unless the number of dimensions is
    % fixed.

    % The dim parameter means we need to work along a specific dimension.
    % The builtin Matlab function handles cases where dim > ndims. We do
    % not do this, as our error check above has trapped these cases and it
    % is not clear that this behaviour is useful. If it turns out one day
    % to be so, the code here will need to be rewritten.
    
    %d = size(A); d(dim) = 1;  % Preallocate an array of the correct size
    %X = quaternion(zeros(d)); % for the result.
        
    error('prod cannot yet handle arrays of dimension greater than 2.')    
end

function Y = product(X)
% Internal function to compute the product of the elements of a vector.
% We need to use subscripted references, but these do not work inside a
% class method (which this is). See the file 'implementation notes.txt',
% item 8.

% TODO Eliminate the recursion here. Iteration would be faster.

if length(X) == 1
    Y = X;
else
    L = length(X); % We compute this because we cannot use 'end' notation.
    Y = subsref(X, substruct('()', {1})) ... % X(1)
        .* ...
        product(subsref(X, substruct('()', {2:L}))); % product(X(2:end))
end

% $Id: prod.m 1004 2017-11-15 17:14:09Z sangwine $
