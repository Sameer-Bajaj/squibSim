function [flutterMode, flutterSpeed] = readFlutterF06(f06File)
    fid = fopen(f06File, 'r');
    if fid == -1, error('Cannot open file: %s', f06File); end
    cleanupObj = onCleanup(@() fclose(fid));

    data = struct();
    currentPoint = -1;
    
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), break; end

        % identify the "point"
        if contains(line, 'POINT')
            tok = regexp(line, 'POINT\s*=\s*(\d+)', 'tokens', 'once');
            if ~isempty(tok)
                currentPoint = str2double(tok{1});
                pKey = sprintf('p%d', currentPoint);
                if ~isfield(data, pKey)
                    data.(pKey).vel = [];
                    data.(pKey).damp = [];
                    data.(pKey).point = currentPoint;
                end
            end
            continue;
        end

        % 2. extract data
        % look for lines that start with at least 3 numbers (KFREQ, 1/K, VELOCITY)
        nums = sscanf(line, '%f');
        if numel(nums) >= 4 && currentPoint > 0
            pKey = sprintf('p%d', currentPoint);
            % Only grab if it looks like a flutter data row (Velocity is index 3, Damping is 4)
            data.(pKey).vel(end+1)  = nums(3);
            data.(pKey).damp(end+1) = nums(4);
        end
    end

    %% Detect Flutter Crossings
    flutter_modes = [];
    flutter_speeds = [];
    pointKeys = fieldnames(data);

    for k = 1:numel(pointKeys)
        pKey = pointKeys{k};
        v = data.(pKey).vel;
        g = data.(pKey).damp;
        
        if numel(v) < 2, continue; end

        for i = 1:numel(g)-1
            % check for any zero crossing (sign change)
            if sign(g(i)) ~= sign(g(i+1)) && g(i+1) ~= 0
                % interpolate flutter velocity:
                v_flutter = v(i) - g(i) * (v(i+1) - v(i)) / (g(i+1) - g(i));
                flutter_modes(end+1)  = data.(pKey).point;
                flutter_speeds(end+1) = v_flutter;
                break; 
            end
        end
    end
    if isempty(flutter_speeds)
        error('readFlutterF06:noFlutterFound', 'No flutter found in %s', f06File);
    end

    [flutterSpeed, idx] = min(flutter_speeds);
    flutterMode = flutter_modes(idx);
end