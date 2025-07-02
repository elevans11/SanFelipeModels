function [d, azimuth] = haversine_distance(R,lat1, lon1, lat2, lon2)
     % Expand scalars to match the size of destination arrays
    if isscalar(lat1)
        lat1 = lat1 * ones(size(lat2));
    end
    if isscalar(lon1)
        lon1 = lon1 * ones(size(lon2));
    end

    % Check size compatibility
    if ~isequal(size(lat1), size(lat2), size(lon1), size(lon2))
        error('Inputs must be scalars or arrays of the same size after broadcasting.');
    end

    % Convert degrees to radians
    lat1_rad = deg2rad(lat1);
    lon1_rad = deg2rad(lon1);
    lat2_rad = deg2rad(lat2);
    lon2_rad = deg2rad(lon2);

    % Differences
    dlat = lat2_rad - lat1_rad;
    dlon = lon2_rad - lon1_rad;

    % Haversine formula
    a = sin(dlat/2).^2 + cos(lat1_rad) .* cos(lat2_rad) .* sin(dlon/2).^2;
    c = 2 .* atan2(sqrt(a), sqrt(1 - a));
    d = R .* c;

    % Azimuth calculation
    y = sin(dlon) .* cos(lat2_rad);
    x = cos(lat1_rad) .* sin(lat2_rad) - ...
        sin(lat1_rad) .* cos(lat2_rad) .* cos(dlon);
    azimuth_rad = atan2(y, x);

    % Convert azimuth to degrees and normalize
    azimuth = mod(rad2deg(azimuth_rad), 360);
end