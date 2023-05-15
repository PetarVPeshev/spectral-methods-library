function norm_magn_matrix = norm_magnitude(magn_matrix, varargin)
%NORM_MAGNITUDE This function normalizes the magnitude to the maximum
%magnitude
%   The function takes the magnitude (scalar, vector, or matrix) and
%   returns the normalized to the maximum magnitude. The output can be in
%   dB or in linear scale. If scale is not specified (by 'dB' or 'linear'),
%   the default is 'dB' scale.
    if isempty(varargin)
        scale_type = 'dB';
    elseif length(varargin) == 1
        scale_type = varargin{1};
        if ~( strcmp(scale_type, 'dB') || strcmp(scale_type, 'linear') )
            error('Invalid argument');
        end
    else
        error('Invalid argument');
    end
    
    max_magn = max( max( abs(magn_matrix) ) );
    if strcmp(scale_type, 'dB')
        norm_magn_matrix = 20 * log10( abs(magn_matrix) ) - ...
            20 * log10(max_magn);
    else
        norm_magn_matrix = abs(magn_matrix) / max_magn;
    end
end