function FT = ft_current(k, varargin)
%FT_CURRENT This function calculates the current's Fourier Transform
%   Detailed explanation goes here
% , k_comp, width, length, antenna, orientation
    antenna = varargin{end - 1};
    orientation = varargin{end};
    
    if strcmp(antenna, 'dipole')
        k_comp = varargin{1};
        width = varargin{2};
        length = varargin{3};
        dielectric_er = varargin{4};

        dielectric_k = k * sqrt(dielectric_er);

        keq = (k + dielectric_k) / 2;

        if strcmp(orientation, 'x')
            T = sinc( k_comp(:, :, 2) * width / (2 * pi) );
            F = 2 * keq * ( cos(k_comp(:, :, 1) * length / 2) - ...
                cos(keq * length / 2) ) ./ ( (keq^2 - ...
                k_comp(:, :, 1).^2) * sin(keq * length / 2) );

            FT = zeros( [size(k_comp, 1, 2), 2] );
            FT(:, :, 1) = F .* T;
        else
            error('Not implemented');
        end
    elseif strcmp(antenna, 'circular')
        a = varargin{1};
        theta = varargin{2};

        J = 2 * pi * (a ^ 2) * besselj(1, k * a * sin(theta)) ./ ...
            (k * a * sin(theta));

        FT = zeros( [size(J, 1, 2), 3] );
        
        if strcmp(orientation, 'y')
            FT(:, :, 2) = J;
        else
            error('Not implemented');
        end
    end
end

