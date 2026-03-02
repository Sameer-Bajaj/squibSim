function runNastran(bdfFile, nastranExe)
    if nargin < 2
        nastranExe = "C:\Program Files\MSC.Software\NaPa_SE\20241\Nastran\bin\nastran.exe";
    end
    % Construct the expected F06 filename
    [fDir, fName, ~] = fileparts(bdfFile);
    f06File = fullfile(fDir, [fName, '.f06']);

    % Delete old F06 if it exists to avoid false positives
    if exist(f06File, 'file'), delete(f06File); end

    % Build the command. 
    % Added 'batch=no' or similar keywords can sometimes help force blocking,
    % but the system() call below is the primary driver.
    cmd = sprintf('"%s" "%s" out="%s" news=no scratch=no', nastranExe, bdfFile, fDir);
   
    fprintf('Running: %s\n', fName);
    
    % Execute Nastran
    [status, ~] = system(cmd);
    
    if status ~= 0
        error('runNastran:failed', 'Nastran launcher failed with status %d.', status);
    end

    %% wait for .f06 to exist and finish writing
    timeout = 180; % 3 minutes max wait
    elapsed = 0;
    pauseTime = 2; % Check every 2 seconds
    
    fprintf('Waiting for %s.f06 to stabilize...\n', fName);
    
    while elapsed < timeout
        if exist(f06File, 'file')
            % Even if it exists, Nastran might still be writing to it.
            % We try to open it with append access. If it fails, it's locked.
            fid = fopen(f06File, 'a');
            if fid ~= -1
                fclose(fid);
                fprintf('Nastran run complete and file released.\n');
                return; % SUCCESS
            end
        end
        
        pause(pauseTime);
        elapsed = elapsed + pauseTime;
    end

    error('runNastran:timeout', 'Timed out waiting for %s to finish.', f06File);
end