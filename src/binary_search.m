function [found, pos] = binary_search(arr, elem)
% binary_search function performs binary search over a sorted array looking
% for a specific object given by elem. It returns whether the element could
% be found or not, and the position where it is or where it should be
%
% INPUT:
%   - arr: array over which we want to perform binary search
%   - elem: element that we are interested in finding
%
% OUTPUT:
%   - found: boolean that tells us whether or not the element could be
%       found in the array
%   - pos: integer that will have the position where the element was found
%       (if it could be found) or the position where it should be to keep
%       the array sorted

    found = false;
    st_pos = 1;
    end_pos = size(arr, 2);
    
    while (st_pos <= end_pos && ~found)
        half = floor((st_pos + end_pos) / 2);
        if elem < arr(half)
            end_pos = half - 1;
        elseif elem > arr(half)
            st_pos = half + 1;
        else
            found = true;
        end
    end
    
    if found
        pos = half;
    else
        pos = st_pos;
    end
end

