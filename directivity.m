function [ dir, rad_intensity, rad_power ] = directivity(relat_permit, ...
    field, sph_grid, r)
%DIRECTIVITY This function calculates the directivity, radiation intensity,
% and radiated power
%   Detailed explanation goes here
    wave_impedance = 376.730313668 / sqrt(relat_permit);

    dth = sph_grid(1, 2, 1) - sph_grid(1, 1, 1);
    dph = sph_grid(2, 1, 2) - sph_grid(1, 1, 2);
    field_total = total_field(field);
    theta = sph_grid(:, :, 1);
    
    rad_intensity = (field_total .^ 2) .* (r .^ 2) / (2 * wave_impedance);
    [row, col] = find( isnan(rad_intensity) );
    if ~isempty(col)
        warning(['Found singularities (NaN) in radiation intensity; ' ...
            'set NaN cells to zero in radiation intensity calculation.']);
        rad_intensity(row, col) = 0;
    end
    rad_power = sum( sum(rad_intensity .* sin(theta)) ) * dth * dph;
    dir = rad_intensity * 4 * pi / rad_power;
end
