function field_total = total_field(field, varargin)
%TOTAL_FIELD This function calculates the total field magnitude from the
%field's components
%   Detailed explanation goes here
    if isempty(varargin)
        field_total = sqrt( abs(field(:, :, 1)) .^ 2 ...
            + abs(field(:, :, 2)) .^ 2 + abs(field(:, :, 3)) .^ 2);
    elseif strcmp(varargin{1}, '5D')
        field_total = sqrt( abs(field(:, :, 1, :, :)) .^ 2 ...
            + abs(field(:, :, 2, :, :)) .^ 2 ...
            + abs(field(:, :, 3, :, :)) .^ 2);
        field_total = permute(field_total, [1 2 4 5 3]);
    else
        error('Error. Not a valid input.');
    end
end

