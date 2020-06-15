classdef HashTable
% Honestly, name is fairly self-explanatory
    
    properties
        size % original size of hashtable
        keys % vector of keys (which are hashed values)
        values % cell array of values 
    end
    
    methods
        function obj = HashTable(order)
            % Construct an instance of this class
            % Generates table of length 2^order
            if nargin == 0
                order = 20;
            end
            obj.size = 2^order;
            obj.keys = zeros(1, obj.size);
            obj.values = cell(1, obj.size);
        end
        
        function obj = putElement(obj, element, group)
            % adds hash of element to the table of keys
            % hash is based on the group to which element belongs
            hash = group.hash(element);
            ind = bitand(hash, obj.size - 1) + 1;
            while obj.keys(ind) ~= 0
                    ind = ind + 1;
                    if ind >= obj.size
                        ind = 1;
                    end
            end
            obj.keys(ind) = hash;
        end
        
        function bool = containsElement(obj, element, group)
            % checks for presence of hash of element in obj.keys
            % returns 1 if hash is present and 0 otherwise
            hash = group.hash(element);
            ind = bitand(hash, obj.size - 1) + 1;
            stop_ind = ind - 1;
            bool = 0;
            while obj.keys(ind) ~= 0
                if obj.keys(ind) == hash
                    bool = 1;
                    break
                elseif ind >= obj.size
                    ind = 1;
                elseif ind == stop_ind
                    err('Error: key not found')
                    break
                else
                    ind = ind + 1;
                end
            end
        end
        
        function obj = put(obj, key, value)
            % adds value to table based on index of key 
            % key must be an object that can be hashed by function Hash
            hash = Hash(key);
            ind = bitand(hash, obj.size - 1) + 1;
            while obj.keys(ind) ~= 0
                    ind = ind + 1;
                    if ind >= obj.size
                        ind = 1;
                    end
            end
            obj.keys(ind) = hash;
            obj.value{ind} = value;
        end
        
                
        function bool = contains(obj, key)
            % checks for presence of key
            hash = Hash(key);
            ind = bitand(hash, obj.size - 1) + 1;
            bool = 0;
            while obj.keys(ind) ~= 0
                if obj.keys(ind) == hash
                    bool = 1;
                    break
                elseif ind >= obj.size
                    ind = 1;
                elseif ind == stop_ind
                    err('Error: key not found')
                    break
                else
                    ind = ind + 1;
                end
            end
        end
        
        function value = get(obj, key)
            % checks for presence of a value corresponding to key
            hash = Hash(key);
            ind = bitand(hash, obj.size - 1) + 1;
            while obj.keys(ind) ~= 0
                if obj.keys(ind) == hash
                    value = obj.values{ind};
                    break
                elseif ind >= obj.size
                    ind = 1;
                elseif ind == stop_ind
                    err('Error: key not found')
                    break
                else
                    ind = ind + 1;
                end
            end
        end
        
        function value = getOrDefault(obj, key, valueIfNotThere)
            % checks for presence of a value corresponding to key
            hash = Hash(key);
            ind = bitand(hash, obj.size - 1) + 1;
            while obj.keys(ind) ~= 0
                if obj.keys(ind) == hash
                    value = obj.values{ind};
                    break
                elseif ind >= obj.size
                    ind = 1;
                elseif ind == stop_ind
                    value = valueIfNotThere;
                    break
                else
                    ind = ind + 1;
                end
            end
        end
        
        function value = getOrElse(obj, key, funcHandle)
            % checks for presence of a value corresponding to key
            hash = Hash(key);
            ind = bitand(hash, obj.size - 1) + 1;
            while obj.keys(ind) ~= 0
                if obj.keys(ind) == hash
                    value = obj.values{ind};
                    break
                elseif ind >= obj.size
                    ind = 1;
                elseif ind == stop_ind
                    value = funcHandle(key);
                    break
                else
                    ind = ind + 1;
                end
            end
        end
        
        
        function value = getIfAbsent(obj, key, howtoCompute)
            % checks for presence of a value corresponding to key
            hash = Hash(key);
            ind = bitand(hash, obj.size - 1) + 1;
            while obj.keys(ind) ~= 0
                if obj.keys(ind) == hash
                    value = obj.values{ind};
                    break
                elseif ind >= obj.size
                    ind = 1;
                elseif ind == stop_ind
                    [new_key, value] = howtoCompute(key);
                    obj.put(new_key, value)
                    break
                else
                    ind = ind + 1;
                end
            end
        end
        
    end
end

