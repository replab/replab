function kv1 = map(kv, funs)
% Maps the values in the given key/value pairs through functions
%
% For each key/value pair where the key is a field of the struct ``funs``, this function
% maps the corresponding value through the function in the struct.
%
% Example:
%   >>> kv = {'a', 1, 'b', 2};
%   >>> funs = struct('a', @(x) 3*x);
%   >>> kv1 = replab.kv.map(kv, funs)
%       {'a', 3, 'b', 2}
% Args:
%   kv (cell(1,\*)): Key/value pairs
%   fun (struct): Struct whose value are function handles
%
% Returns:
%   cell(1,\*): The updated key/value pairs
    kv1 = kv;
    for i = 1:2:length(kv)
        if isfield(funs, kv{i})
            fun = funs.(kv{i});
            kv1{i+1} = fun(kv{i+1});
        end
    end
end
