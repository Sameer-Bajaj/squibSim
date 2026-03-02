%{

SquibSim™

By: Sameer Bajaj

Allows for rapid solution of a large set of flutter problems.

Directory "mainLoc" should include folder(s) containing exported BDFs from 
ANSYS Modal, with contacts having MPC formulation and node matching between
split fin halves.

%}

clear all; close all; clc;

%% BDF Write Input Parameters
mach = 2.0;
relDensity = 0.6;
LMODES = 6;
numSpanGP = 8;
numChordGP = 16;
NBOX = 11; % for machBox
vHigh = 5000;
vLow = 1200;
mainLoc = 'C:\Users\Owner\Downloads\squibSimDir';
numInstances = 2;
numSteps = 47; % number of velocity steps
strucDampBool = true;
SDAMP = 0.005;
redFreq = [0.01, 0.02, 0.05, 0.10, 0.2, 0.25]; % list of reduced freq

% Defaults, DO NOT COMMENT OUT!::
defReqFreq = [0.01, 0.02, 0.05, 0.10, 0.2, 0.25];
% =========================================================================

% This next segment handles the user input.
%
% You can choose to either directly input the variables above, or comment
% those lines out and be prompted to enter your parameters each time.

% =========================================================================

if ~exist('mainLoc', 'var')
    mainLoc = input(sprintf('%-66s', 'Enter your directory [C:\Users\...]:'), 's');
end

if ~exist('numInstances', 'var')
    numInstances = input(sprintf('%-64s', 'Enter the maximum number of NASTRAN instances:'));
end

if ~exist('mach', 'var')
    mach = input(sprintf('%-64s', 'Enter your free stream mach (M > 1):'));
    while mach <= 1
        warning('Flow must be supersonic!')
        mach = input(sprintf('%-64s', 'Enter your free stream mach (M > 1):'));
    end
end

if ~exist('relDensity', 'var')
    relDensity = input(sprintf('%-64s', 'Enter relative density (rho/rho_0):'));
end

if ~exist('vLow', 'var') || ~exist('vHigh', 'var')
    vLow = input(sprintf('%-64s', 'Enter your low bound velocity (ft/s):'));
    vHigh = input(sprintf('%-64s', 'Enter your high bound velocity (ft/s):'));
end

if ~exist('numSteps', 'var')
    numSteps = input(sprintf('%-64s', 'Enter the number of velocity steps:'));
end

if ~exist('LMODES', 'var')
    LMODES = input(sprintf('%-64s', 'Enter number of modes to include (LMODES):'));
end

if ~exist('numSpanGP', 'var') || ~exist('numChordGP', 'var')
    numSpanGP = input(sprintf('%-64s', 'Enter number of spanwise grid points:'));
    numChordGP = input(sprintf('%-64s', 'Enter number of chordwise grid points:'));
end

if ~exist('NBOX', 'var') && mach <= 2.4
    NBOX = input(sprintf('%-64s', 'Enter number of Mach Boxes:'));
end
if ~exist('strucDampBool', 'var')
    resp = input(sprintf('%-64s', 'Include structural damping? (y/n):'), 's');
    strucDampBool = strcmpi(resp, 'y');
end

if strucDampBool && ~exist('SDAMP', 'var')
    SDAMP = input(sprintf('%-64s', 'Enter structural damping value (G):'));
end

if ~exist('redFreq', 'var')
    disp('User did not specify reduced frequencies!')
    fprintf(['Using default reduced frequency values: [%.3f, %.3f, ' ...
        '%.3f, %.3f, %.3f, %.3f]\n'], defReqFreq(1), defReqFreq(2), ...
        defReqFreq(3), defReqFreq(4), defReqFreq(5), defReqFreq(6))
    redFreq = defReqFreq;
end

files = dir(mainLoc);

% exclude non-folders:
foldFlag = [files.isdir];
files = files(foldFlag);
fileNames = {files.name};
folders = files(~((strcmp(fileNames, '.') | strcmp(fileNames, '..'))));

foldNames = {folders.name}; % gets names of all folders in dir
numTrials = numel(foldNames);
numGroups = ceil((numTrials - mod(numTrials, numInstances)) / numInstances) ...
    + (mod(numTrials, numInstances) ~= 0); % number of groups to run

%% Parallel Execution and Parsing

% 1. setup parallel pool
p = gcp('nocreate');
if isempty(p)
    parpool(numInstances);
elseif p.NumWorkers ~= numInstances
    delete(p);
    parpool(numInstances);
end

% 2. initialize results containers
allFlutterSpeeds = NaN(numTrials, 1);
allFlutterModes  = NaN(numTrials, 1);

parfor folderIDX = 1:numTrials
    % Staggered start (first batch only) ---
    if folderIDX <= numInstances
        pause(4 * (folderIDX - 1)); 
    end
    
    foldLoc = fullfile(mainLoc, foldNames{folderIDX});
    maxRetries = 12;
    success = false;

    % Retry until solution converges
    for attempt = 1:maxRetries
        try
            % a) write the BDF
            if mach <= 2.4
                bdfFile = writeMachBoxBDF(mach, foldLoc, vLow, vHigh, numSteps, LMODES, ...
                    SDAMP, strucDampBool, relDensity, numSpanGP, numChordGP, redFreq, NBOX);
            else
                bdfFile = writePistonBDF(mach, foldLoc, vLow, vHigh, numSteps, LMODES, ...
                    SDAMP, strucDampBool, relDensity, numSpanGP, numChordGP, redFreq);
            end

            % b) run NASTRAN
            runNastran(bdfFile); 

            % c: parse  the .f06 file
            [fPath, fName, ~] = fileparts(bdfFile);
            f06File = fullfile(fPath, [fName, '.f06']);
            
            if exist(f06File, 'file')
                % d. get flutter mode and velocity
                [m, v] = readFlutterF06(f06File);
                allFlutterModes(folderIDX)  = m;
                allFlutterSpeeds(folderIDX) = v;
                
                fprintf('Trial [%s]: Success on attempt %d. Speed %.2f ft/s\n', ...
                    foldNames{folderIDX}, attempt, v);
                
                success = true; 
                break; % EXIT the retry loop
            else
                error('SquibSim:f06Missing', 'F06 file not found after run.');
            end
            
        catch ME
            fprintf('Trial [%s]: Attempt %d failed. Error: %s\n', ...
                foldNames{folderIDX}, attempt, ME.message);
                
            if attempt < maxRetries
                pause(5 * attempt); % wait exponentially longer each retry
            end
        end
    end
    
    if ~success
        fprintf('Trial [%s]: Critical Failure. All %d attempts exhausted.\n', ...
            foldNames{folderIDX}, maxRetries);
    end
end

%% Final Summary Table
% creates header (bunch of repeated  (repmat) equal signs)
fprintf(['\n', repmat('=', 1, 50), '\n']); 
fprintf('%-20s | %-10s | %-12s\n', 'Trial Name', 'Mode', 'Speed (ft/s)');
fprintf([repmat('-', 1, 50), '\n']);

for i = 1:numTrials
    if ~isnan(allFlutterSpeeds(i))
        fprintf('%-20s | %-10d | %-12.2f\n', foldNames{i}, allFlutterModes(i), allFlutterSpeeds(i));
    else
        fprintf('%-20s | %-10s | %-12s\n', foldNames{i}, 'N/A', 'Stable');
    end
end

%% Post-Run Root Directory Cleanup
fprintf('\nCleaning up root directory... ');

% get the path where this script is actually located
curLoc = pwd; % current program location
rootFiles = dir(curLoc);
% delete everything that's not a matlab file:
for i = 1:numel(rootFiles)
    [~, fName, fExt] = fileparts(rootFiles(i).name);
    % skip directories
    if rootFiles(i).isdir, continue; end
    % define safe filetypes
    isMatlabFile = strcmpi(fExt, '.m') || strcmpi(fExt, '.mat') || strcmpi(fExt, '.p');
    
    if ~isMatlabFile
        try
            delete(fullfile(curLoc, rootFiles(i).name));
        catch
        end
    end
end

fprintf('Done.\n');
