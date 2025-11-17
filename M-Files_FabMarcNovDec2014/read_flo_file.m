function [u, v] = read_flo_file(filename)
    % READ_FLO_FILE Read a .flo optical flow file (Middlebury format)
    % Usage: [u, v] = read_flo_file('flow.flo')

    fid = fopen(filename, 'rb');
    if fid < 0
        error('Could not open %s', filename);
    end

    tag = fread(fid, 1, 'float32');
    if tag ~= 202021.25
        fclose(fid);
        error('Invalid .flo file (wrong tag: %f)', tag);
    end

    width  = fread(fid, 1, 'int32');
    height = fread(fid, 1, 'int32');

    % Read flow data (interleaved u and v)
    data = fread(fid, [2, width * height], 'float32');
    fclose(fid);

    % Reshape and separate u and v
    data = reshape(data, [2, width, height]);
    data = permute(data, [3 2 1]);  % (height, width, 2)
    u = data(:, :, 1);
    v = data(:, :, 2);
end