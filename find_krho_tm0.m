function krho_tm0 = find_krho_tm0(k0, varargin)
%FIND_KRHO_TM0 Summary of this function goes here
%   Detailed explanation goes here
    if strcmp(varargin{1}, 'SemiInfiniteSuperstrate')
        air_length = varargin{2};
        slab_er = varargin{3};
    
        % Approximated kz
        air_kz_tm0 = sqrt(1j * k0 / (sqrt(slab_er) * air_length));
    
        % Approximated krho in air
        krho_tm0 = sqrt(k0 ^ 2 - air_kz_tm0 ^ 2);
            
        % Step k
        delta_k = k0 / 500;
            
        % Newton's method for TM krho
        krho_tm0_prev = 0;
        idx = 0;
        while abs(krho_tm0 - krho_tm0_prev) > 0.00001 && idx < 500
            [~, D_tm] = dispersion_eqn(k0, krho_tm0, ...
                'SemiInfiniteSuperstrate', air_length, slab_er);
        
            [~, D_tm_1] = dispersion_eqn(k0, krho_tm0 + delta_k / 2, ...
                'SemiInfiniteSuperstrate', air_length, slab_er);
            [~, D_tm_2] = dispersion_eqn(k0, krho_tm0 - delta_k / 2, ...
                'SemiInfiniteSuperstrate', air_length, slab_er);
        
            D_tm_prime = (D_tm_1 - D_tm_2) / delta_k;
            krho_tm0_prev = krho_tm0;
            krho_tm0 = krho_tm0_prev - D_tm / D_tm_prime;
            idx = idx + 1;
        end
    end
end

