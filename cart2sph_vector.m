function sph_vector = cart2sph_vector(cart_vector, varargin)
%CART2SPH_VECTOR This function converts spherical vector in cartesian
%coordinate system to a vector in spherical coordinate system
%   The function takes M-by-N-by-3 matrix for the cartesian vector, and 
%   M-by-N matrix for Theta and Phi. The cartesian matrix's third dimension
%   represents the cartesian unit vectors as x, y, and z respectively. The 
%   M-by-N matricies of each cartesian unit vector holds its magnitude for 
%   the corresponding vector location at the elevation and azimuth. The 
%   function returns M-by-N-by-3 matrix for the spherical vector with third
%   dimensions representing the radial, elevation angle, and azimuth angle 
%   unit vectors respectively.
    if length(varargin) == 4
        for idx = 1:2:4
            if strcmp(varargin{idx}, 'Theta')
                theta = varargin{idx + 1};
            else
                phi = varargin{idx + 1};
            end
        end
    elseif length(varargin) == 2
        theta = varargin{1};
        phi = varargin{2};
    elseif length(varargin) == 1
        theta = varargin{1}(:, :, 1);
        phi = varargin{1}(:, :, 2);
    else
        error('Invalid arguments');
    end

    % Spherical vector computation
    sph_vector = NaN( [size(cart_vector)] );
    sph_vector(:, :, 1) = ...
        cart_vector(:, :, 1) .* sin(theta) .* cos(phi) + ...
        cart_vector(:, :, 2) .* sin(theta) .* sin(phi) + ...
        cart_vector(:, :, 3) .* cos(theta);
    sph_vector(:, :, 2) = ...
        cart_vector(:, :, 1) .* cos(theta) .* cos(phi) + ...
        cart_vector(:, :, 2) .* cos(theta) .* sin(phi) - ...
        cart_vector(:, :, 3) .* sin(theta);
    sph_vector(:, :, 3) = ...
        - cart_vector(:, :, 1) .* sin(phi) + ...
        cart_vector(:, :, 2) .* cos(phi);
end