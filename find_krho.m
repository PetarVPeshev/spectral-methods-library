function [krho_te, krho_tm] = find_krho(k0, krho, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if strcmp(varargin{1}, 'GroundSlab')
        slab_length = varargin{2};
        dielectric_er = varargin{3};

        [D_te, D_tm] = dispersion_eqn(k0, krho, ...
            'GroundSlab', slab_length, dielectric_er);
    
        % Peaks
        [~, peak_te] = findpeaks( abs(1 ./ D_te) );
        [~, peak_tm] = findpeaks( abs(1 ./ D_tm) );
        
        % Guess krho
        krho_te = krho(peak_te);
        krho_tm = krho(peak_tm);
        
        % Step k
        delta_k = k0 / 500;
        
        % Newton's method for TE krho
        krho_te_prev = 0;
        if isempty(krho_te)
            krho_te = k0;
            krho_te_prev = k0;
        end

        if length(krho_te) ~= 1
            krho_te = krho_te(end);
            warning('Several solutions to dispersion equation are found.');
        end

        while abs(krho_te - krho_te_prev) > 0.00001
            [D_te, ~] = dispersion_eqn(k0, krho_te, ...
                'GroundSlab', slab_length, dielectric_er);
        
            [D_te_1, ~] = dispersion_eqn(k0, krho_te + delta_k / 2, ...
                'GroundSlab', slab_length, dielectric_er);
            [D_te_2, ~] = dispersion_eqn(k0, krho_te - delta_k / 2, ...
                'GroundSlab', slab_length, dielectric_er);
        
            D_te_prime = (D_te_1 - D_te_2) / delta_k;
            krho_te_prev = krho_te;
            krho_te = krho_te - D_te / D_te_prime;
            
            if imag(krho_te) ~= 0
                krho_te = k0;
                break;
            end
        end
        
        % Newton's method for TM krho
        krho_tm_prev = 0;
        if isempty(krho_tm)
            krho_tm = k0;
            krho_tm_prev = k0;
        end

        if length(krho_tm) ~= 1
            krho_tm = krho_tm(end);
            warning('Several solutions to dispersion equation are found.');
        end
        
%         while abs(krho_tm - krho_tm_prev) > 0.00001 && ~isempty(krho_tm)
        while abs(krho_tm - krho_tm_prev) > 0.00001
            [~, D_tm] = dispersion_eqn(k0, krho_tm, ...
                'GroundSlab', slab_length, dielectric_er);
        
            [~, D_tm_1] = dispersion_eqn(k0, krho_tm + delta_k / 2, ...
                'GroundSlab', slab_length, dielectric_er);
            [~, D_tm_2] = dispersion_eqn(k0, krho_tm - delta_k / 2, ...
                'GroundSlab', slab_length, dielectric_er);
        
            D_tm_prime = (D_tm_1 - D_tm_2) / delta_k;
            krho_tm_prev = krho_tm;
            krho_tm = krho_tm - D_tm / D_tm_prime;
            
            if imag(krho_tm) ~= 0
                krho_tm = k0;
                break;
            end
        end
    elseif strcmp(varargin{1}, 'Superstrate')
        air_length = varargin{2};
        slab_length = varargin{3};
        slab_er = varargin{4};
        
        wavelength = 2 * pi / k0;
        h_bar = air_length / wavelength;
    
        % Approximated kz
        air_kz_te = k0 * (2 * pi * h_bar * slab_er + 1j) ...
            / ((pi * slab_er * (2 * h_bar) ^ 2) ...
            * (1 + 1 / ( (2 * pi * h_bar * slab_er) ^ 2)));
        air_kz_tm = (k0 / (4 * h_bar)) ...
            * (1 + sqrt(1 + 8j * h_bar / (pi * slab_er)));
    
        % Approximated krho
        krho_te = sqrt(k0 ^ 2 - air_kz_te ^ 2);
        krho_tm = sqrt(k0 ^ 2 - air_kz_tm ^ 2);
            
        % Step k
        delta_k = k0 / 500;
            
        % Newton's method for TE krho
        krho_te_prev = 0;
        idx = 0;
        while abs(krho_te - krho_te_prev) > 0.00001 && idx < 500
            [D_te, ~] = dispersion_eqn(k0, krho_te, ...
                'Superstrate', air_length, slab_length, slab_er);
        
            [D_te_1, ~] = dispersion_eqn(k0, krho_te + delta_k / 2, ...
                'Superstrate', air_length, slab_length, slab_er);
            [D_te_2, ~] = dispersion_eqn(k0, krho_te - delta_k / 2, ...
                'Superstrate', air_length, slab_length, slab_er);
        
            D_te_prime = (D_te_1 - D_te_2) / delta_k;
            krho_te_prev = krho_te;
            krho_te = krho_te_prev - D_te / D_te_prime;
            idx = idx + 1;
        end
            
        % Newton's method for TM krho
        krho_tm_prev = 0;
        idx = 0;
        while abs(krho_tm - krho_tm_prev) > 0.00001 && idx < 500
            [~, D_tm] = dispersion_eqn(k0, krho_tm, ...
                'Superstrate', air_length, slab_length, slab_er);
        
            [~, D_tm_1] = dispersion_eqn(k0, krho_tm + delta_k / 2, ...
                'Superstrate', air_length, slab_length, slab_er);
            [~, D_tm_2] = dispersion_eqn(k0, krho_tm - delta_k / 2, ...
                'Superstrate', air_length, slab_length, slab_er);
        
            D_tm_prime = (D_tm_1 - D_tm_2) / delta_k;
            krho_tm_prev = krho_tm;
            krho_tm = krho_tm_prev - D_tm / D_tm_prime;
            idx = idx + 1;
        end
    elseif strcmp(varargin{1}, 'SemiInfiniteSuperstrate')
        air_length = varargin{2};
        slab_er = varargin{3};
        
        wavelength = 2 * pi / k0;
        h_bar = air_length / wavelength;
    
        % Approximated kz
        air_kz_te = k0 * (2 * pi * h_bar * sqrt(slab_er) + 1j) ...
            / ((pi * sqrt(slab_er) * (2 * h_bar) ^ 2) ...
            * (1 + 1 / ( slab_er * (2 * pi * h_bar) ^ 2)));
        air_kz_tm = (k0 / (4 * h_bar)) ...
            * (1 + sqrt(1 + 8j * h_bar / (pi * sqrt(slab_er))));
    
        % Approximated krho in air
        krho_te = sqrt(k0 ^ 2 - air_kz_te ^ 2);
        krho_tm = sqrt(k0 ^ 2 - air_kz_tm ^ 2);
            
        % Step k
        delta_k = k0 / 500;
            
        % Newton's method for TE krho
        krho_te_prev = 0;
        idx = 0;
        while abs(krho_te - krho_te_prev) > 0.00001 && idx < 500
            [D_te, ~] = dispersion_eqn(k0, krho_te, ...
                'SemiInfiniteSuperstrate', air_length, slab_er);
        
            [D_te_1, ~] = dispersion_eqn(k0, krho_te + delta_k / 2, ...
                'SemiInfiniteSuperstrate', air_length, slab_er);
            [D_te_2, ~] = dispersion_eqn(k0, krho_te - delta_k / 2, ...
                'SemiInfiniteSuperstrate', air_length, slab_er);
        
            D_te_prime = (D_te_1 - D_te_2) / delta_k;
            krho_te_prev = krho_te;
            krho_te = krho_te_prev - D_te / D_te_prime;
            idx = idx + 1;
        end
            
        % Newton's method for TM krho
        krho_tm_prev = 0;
        idx = 0;
        while abs(krho_tm - krho_tm_prev) > 0.00001 && idx < 500
            [~, D_tm] = dispersion_eqn(k0, krho_tm, ...
                'SemiInfiniteSuperstrate', air_length, slab_er);
        
            [~, D_tm_1] = dispersion_eqn(k0, krho_tm + delta_k / 2, ...
                'SemiInfiniteSuperstrate', air_length, slab_er);
            [~, D_tm_2] = dispersion_eqn(k0, krho_tm - delta_k / 2, ...
                'SemiInfiniteSuperstrate', air_length, slab_er);
        
            D_tm_prime = (D_tm_1 - D_tm_2) / delta_k;
            krho_tm_prev = krho_tm;
            krho_tm = krho_tm_prev - D_tm / D_tm_prime;
            idx = idx + 1;
        end
    end
end