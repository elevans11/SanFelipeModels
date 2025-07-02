function [lat2, lon2] = reckon_sphere(R,lat1, lon1, azimuth, distance)

    % Convert inputs to radians
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    azimuth = deg2rad(azimuth);

    % Angular distance
    delta = distance / R;

    % Calculate new latitude
    lat2 = asin(sin(lat1) * cos(delta) + cos(lat1) * sin(delta) * cos(azimuth));

    % Calculate new longitude
    lon2 = lon1 + atan2(sin(azimuth) * sin(delta) * cos(lat1), ...
                        cos(delta) - sin(lat1) * sin(lat2));

    % Convert outputs to degrees
    lat2 = rad2deg(lat2);
    lon2 = rad2deg(lon2);

    % Normalize longitude to [-180, 180]
    lon2 = mod(lon2 + 180, 360) - 180;
end
