function field = farfield(k, r, sph_grid, kz, z, sgf, current_ft, varargin)
%FARFIELD This function evaluates the electric or magnetic far-field in 
%spherical coordinates for a given Fourier Transform of electric or 
%magnetic current density
%   Detailed explanation goes here
    z_source = 0;
    if ~isempty(varargin)
        z_source = varargin{1};
    end
    % both sgf and current_ft must be 2D
    sgf_current = permute( pagemtimes( permute(sgf, [3, 4, 1, 2]), ...
        permute(current_ft, [3, 4, 1, 2]) ), [3, 4, 1, 2]);
    % z is the observation point, in reallity it should be 
    % exp(1j * kz .* abs(z - z_source)), where z_source is the location of 
    % the source on the z-axis
    % TODO: add additional varargin to check if the user wants to change
    % the z-axis location of the source
    field = 1j * cart2sph_vector(kz .* sgf_current, sph_grid) ...
        .* exp(1j * kz .* abs(z - z_source)) .* exp(-1j * k * r) ...
        ./ (2 * pi * r);
end

