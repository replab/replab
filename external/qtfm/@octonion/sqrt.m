function Y = sqrt(X)
% SQRT   Square root.
% (Octonion overloading of standard Matlab function.)

% Copyright (c) 2013 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

narginchk(1, 1), nargoutchk(0, 1)

if isreal(X)
    
    % X is a real octonion, and we compute the square root of an
    % isomorphic complex number using the standard Matlab square root
    % function, then construct an octonion with the same axis as the
    % original.
    
    Y = isooctonion(sqrt(isocomplex(X)), X);
else
    
    % X is a complex octonion, and therefore we cannot use the method
    % above for real octonions, because it is not possible to construct
    % an isomorphic complex number. Therefore we use polar form and halve
    % the argument. Note that the modulus and argument here are complex,
    % so the square root of the modulus is complex.
    
    % Preconstruct a result array of zero octonions, of the same size as
    % X and with components of the same class.
    
    Y = zeroso(size(X), class(part(X, 2)));
    
    % If the vector part of any element of X is zero, computing its axis
    % will result in an undefined axis warning. We avoid this warning by
    % computing just the square root of the scalar part in these cases.
    % There may be no such cases, or all the elements in X may have
    % undefined axis.
    
    undefined = abs(normo(v(X))) < eps;
    defined   = ~undefined;
        
    % In order to perform the subscripted assignment using the logical 
    % array undefined we have to use subsasgn and substruct because normal
    % indexing does not work here in a class method.
    % See the file 'implementation notes.txt', item 8, for more details.

    if nnz(undefined) > 0
        % There are some cases of undefined axis.
        U = scalar(subsref(X, substruct('()',  {undefined})));
        Y = subsasgn(Y, substruct('()', {undefined}), octonion(sqrt(U)));
    end
    
    if nnz(defined) > 0
        % There are some cases where the axis is defined.
        D = subsref(X, substruct('()', {defined}));    
        Y = subsasgn(Y, substruct('()', {defined}), ...
                    sqrt(abs(D)) .* exp(axis(D) .* angle(D) ./ 2));
    end
end;

% TODO Consider whether a better algorithm could be devised for the complex
% case based on the Cartesian form. See the quaternion version of this
% function for the full note.

% $Id: sqrt.m 1004 2017-11-15 17:14:09Z sangwine $
