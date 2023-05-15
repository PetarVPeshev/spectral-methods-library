function grid = meshgrid_comb(varargin)
%MESHGRID_COMB Wrap function of meshgrid, which returns the grid matricies
%together in a 3D or 4D matrix
%   The function takes either 2 or 3 vectors. For 2 vectors, it returns
%   3D matrix with the third dimension holding the vector grids in the
%   order of parsing them to the function. For 3 vectors, it returns 4D
%   matrix with the forth dimension holding the vector grids in the order
%   of parsing them to the function.
    x = varargin{1};
    y = varargin{2};
    if length(varargin) == 2
        grid = NaN(length(y), length(x), 2);
        [grid(:, :, 1), grid(:, :, 2)] = meshgrid(x, y);
    elseif length(varargin) == 3
        z = varargin{3};
        grid = NaN(length(y), length(x), length(z), 3);
        [grid(:, :, :, 1), grid(:, :, :, 2), ...
            grid(:, :, :, 3)] = meshgrid(x, y, z);
    else
        error('Invalid arguments');
    end
end

