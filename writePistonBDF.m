%{

NASTRAN BDF WRITER: PISTON METHOD

By: Sameer Bajaj

Function for use by SquibSim™

Intended for use with mach number > 2.4

%}


function bdfFile = writePistonBDF(mach, fileLoc, velLow, velHigh, numSteps, LMODES ...
    , SDAMP, strucDampBool, relDensity, numSpanGP, numChordGP, redFreq)

velocities = linspace(velLow, velHigh, numSteps); % velocity values in ft/s
numStrips = 8;
ft2in = 12;

%% --- READ BDF GRID DATA ---
filenames = {dir(fileLoc).name};
filenames = lower(filenames(:));
bdfFilename = filenames(endsWith(filenames, '.bdf'));
if isempty(bdfFilename)
    error('writePistonBDF:fileNotFound', ...
          'Export BDF from Ansys!');
end
bdfFile = sprintf('%s\\%s', foldLoc, bdfFilename{1});
name = bdfFilename{1};
name = name(1:end - 4);

fid = fopen(bdfFile, 'r');
data = textscan(fid, '%s', 'Delimiter', '\n'); % Read line by line
fclose(fid);
data = data{1}; 
SOL = str2num(data{3}(end-2:end));
formatted = SOL == 145; % check if file has already been formatted
oldDat = data;

%% --- MAP NODES TO AERO PANELS ---
flag = false;
CHEXFLAG = false;
for i = 1:length(data)
    curCell = data{i};
    if length(curCell) <= 4
        continue
    end
    if isstrprop(curCell(1), 'digit')
        if startsWith(data{i - 1}, 'MPC') || flag
            data{i} = sprintf('%-8s%-8s%s', ' ', ' ', curCell);
            flag = true;
            continue
        elseif startsWith(data{i - 1}, 'CTETRA')
            data{i} = sprintf('%-8s%s', ' ', curCell);
            continue
        elseif startsWith(data{i - 1}, 'CHEXA') || CHEXFLAG
            data{i} = sprintf('%-8s%s', ' ', curCell);
            CHEXFLAG = true;
            continue
        end
    end
    if startsWith(curCell, 'MPC')
        flag = false;
        continue
    elseif startsWith(curCell, 'CHEXA')
        CHEXFLAG = false;
    end
    if strcmp(curCell(1:4), 'GRID')
        X(i) = str2double(curCell(25:32));
        Y(i) = str2double(curCell(33:40));
        Z(i) = str2double(curCell(41:48));
        nodeID(i) = str2double(erase(curCell(9:16), " "));
    end
end

tol = 0.0005;
xLE = min(X);
xTE = max(X);
rootChord = xTE - xLE;
yMin = min(Y);
yMax = max(Y);
zMax = max(Z);
span = yMax - yMin;
% find the wedge angle (alpha)
[xWedgeMin, ~] = min(X(Z >= zMax - tol));
[xWedgeMax, ~] = max(X(Z >= zMax - tol));
firWedgeLen = xWedgeMin;
secWedgeLen = xTE - xWedgeMax;
% find tip chord correctly
tipNodes = Y > yMax - 0.005;
if any(tipNodes)
    tipChord = xTE - min(X(tipNodes));
else
    tipChord = rootChord; % fallback if no nodes near tip
end
valid = ~isnan(nodeID) & abs(Z) < tol;
X = X(valid);
Y = Y(valid);
Z = Z(valid);
nodeID = nodeID(valid);

% Y locations of edges of each strip along span
stripEdges = linspace(yMin, yMax, numSpanGP + 1);
stripLen = span / numSpanGP;

% local chord at each strip (linear taper)
chordLensEdge = rootChord - (rootChord - tipChord) * ((stripEdges - yMin) / span);

offset = 0.25; % offset in inches
SET100_nodes = [];

for i = 1:numSpanGP
    lowBoundI = stripEdges(i);
    highBoundI = stripEdges(i + 1);
    inStrip = Y >= lowBoundI & Y < highBoundI;
    
    if ~any(inStrip)
        warning('No nodes in strip %d', i);
        continue
    end
    
    Xs = X(inStrip);
    Ys = Y(inStrip);
    Zs = Z(inStrip);
    nodeIDs = nodeID(inStrip);
    
    xLEs = xTE - chordLensEdge(i);
    
    for j = 1:numChordGP
        lowBoundJ = xLEs + chordLensEdge(i) * (j - 1) / numChordGP;
        highBoundJ = lowBoundJ + chordLensEdge(i) / numChordGP;
        
        chordMidPoint = (lowBoundJ + highBoundJ) / 2;
        
        % find node closest to chord midpoint
        inChordSection = Xs >= lowBoundJ & Xs < highBoundJ;
        if any(inChordSection)
            [~, idx] = min(abs(Xs(inChordSection) - chordMidPoint));
            candidateIDs = nodeIDs(inChordSection);
            ID = candidateIDs(idx);
        else
            % fallback: find nearest node to midpoint
            [~, idx] = min(abs(Xs - chordMidPoint));
            ID = nodeIDs(idx);
        end
        
        % add node if it's not a duplicate
        if isempty(SET100_nodes) || SET100_nodes(end) ~= ID
            SET100_nodes(end + 1) = ID;
        end
    end
end
% Compute average chord len (reference)
refChord = (tipChord + rootChord) / 2;
% Remove any NaN values
SET100_nodes = SET100_nodes(~isnan(SET100_nodes));
% overwrite SET1 header
SET100TAB = {};
SET100TAB{end + 1} = sprintf('%-8s%-8d', 'SET1', 100);
k = 1;
for i = 1:length(SET100_nodes)
    SET100TAB{k} = [SET100TAB{k} sprintf('%-8d', SET100_nodes(i))];

    % line break every 9 entries, tabbed continuation
    if mod(i + 1, 8) == 0
        k = k + 1;
        SET100TAB{end + 1} = '';
        SET100TAB{k} = [SET100TAB{k}, sprintf('%-8s', ' ')];
    end
end
%% --- GET CAERO CARDS ---

sweep = rootChord - tipChord;

% first CAERO card:
x1 = 0; y1 = 0; z1 = 0; x12 = firWedgeLen; x43 = firWedgeLen; x4 = sweep;
y4 = span; z4 = 0;
CAERO101{1} = sprintf('%-8s%-8d%-8d%-8d%-8d%-8s%-8d%-8d', 'CAERO5', 101, ...
    1, 0, numStrips, ' ', 2, 700);
CAERO101{2} = sprintf('%-8s%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f', ...
    ' ', x1, y1, z1, x12, x4, y4, z4, x43);

% second CAERO card:
tipDist = abs(tipChord - firWedgeLen - secWedgeLen); % dist btwn wedges
rootDist = (rootChord - firWedgeLen - secWedgeLen);
x1 = xWedgeMin; y1 = 0; z1 = 0; x12 = rootDist ; x43 = tipDist; x4 = sweep + firWedgeLen;
y4 = span; z4 = 0;
CAERO102{1} = sprintf('%-8s%-8d%-8d%-8d%-8d%-8s%-8d%-8d', 'CAERO5', 201, ...
    1, 0, numStrips, ' ', 2, 701);
CAERO102{2} = sprintf('%-8s%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f', ...
    ' ', x1, y1, z1, x12, x4, y4, z4, x43);

% third CAERO card:
x1 = rootChord - secWedgeLen; y1 = 0; z1 = 0; x12 = secWedgeLen ; x43 = secWedgeLen; 
x4 = x1; y4 = span; z4 = 0;
CAERO103{1} = sprintf('%-8s%-8d%-8d%-8d%-8d%-8s%-8d%-8d', 'CAERO5', 301, ...
    1, 0, numStrips, ' ', 2, 702);
CAERO103{2} = sprintf('%-8s%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f', ...
    ' ', x1, y1, z1, x12, x4, y4, z4, x43);

%% --- COMPUTE THICKNESS INTEGRALS ---
m = zMax / firWedgeLen; % slope of segment 1
NTHICK1 = [m, 0.5*m, 1/3 * m, m^2, 1/2 * m^2, 1/3 * m^2];
NTHICK2str = sprintf('%-8s%-8s%-8s%-8s%-8s%-8s', ...
    '0.', '0.', '0.', '0.', '0.', '0.');
m = zMax / secWedgeLen; % slope of segment 3
NTHICK3 = [-m, -0.5*m, -1/3 * m, m^2, 1/2 * m^2, 1/3 * m^2];
NTHICK3str = sprintf('%-8.4f%-8.4f%-8.4f%-8.4f%-8.4f%-8.4f', ...
    NTHICK3(1), NTHICK3(2), NTHICK3(3), NTHICK3(4), NTHICK3(5), NTHICK3(6));
NTHICK1str = sprintf('%-8.4f%-8.4f%-8.4f%-8.4f%-8.4f%-8.4f', ...
    NTHICK1(1), NTHICK1(2), NTHICK1(3), NTHICK1(4), NTHICK1(5), NTHICK1(6));

%% --- REWRITE HEADER ---
if strucDampBool
    dampLine = '  SDAMP = 2000';
else
    dampLine = ' ';
end
newHeader = {'SOL 145'
    'CEND'
    sprintf('  TITLE = %s', name)
    '  ECHO = NONE'
    dampLine
    '  VECTOR(PLOT)=ALL'
    '  METHOD = 10'
    '  FMETHOD = 30'
    '  MEFFMASS(ALL)=Yes'
    '  SPC = 1'
    'BEGIN BULK'};
SOLIDX = find(startsWith(data, 'SOL'));
bulkIdx = find(startsWith(data, 'BEGIN'));
endIDX = find(startsWith(data, 'END'));
comIDX = find(contains(data, '$ ='), 1, 'first'); % find first comment idx
data(SOLIDX:bulkIdx) = []; % deletes old header
data = [data(1:SOLIDX - 1)
    newHeader
    data(SOLIDX:end)];
%% --- WRITE AEROELASTIC DATA ---

% separation header:
sep = '$ ===============================================================';
% write data that is unchanging
AE_data = {};
if strucDampBool
AE_data{end + 1} = sep;
AE_data{end + 1} = '$ STRUCT DAMPING';
AE_data{end + 1} = sep;
AE_data{end + 1} = ' ';
AE_data{end + 1} = sprintf('%-8s%-8s%-8s', 'PARAM', 'KDAMP', "+1");
AE_data{end + 1} = sprintf('%-8s%-8d', 'TABDMP1', 2000);
AE_data{end + 1} = sprintf('%-8s%-8.1f%-8.3f%-8.1f%-8.3f%-8s', ...
    ' ', 0.0, SDAMP, 1000.0, SDAMP, 'ENDT');
AE_data{end + 1} = ' ';
end
AE_data{end + 1} = sep;
AE_data{end + 1} = '$ AEROELASTIC MODEL PARAMETERS';
AE_data{end + 1} = sep;
AE_data{end + 1} = sprintf('%-8s%-8s%-8s%-8s%-8s', ...
 '$ AERO:', 'ACSID', 'VELOC', 'RefChrd', 'Density');
AE_data{end + 1} = sprintf('%-8s%-8d%-8s%-8.3f%-8s', 'AERO', 0, ' ', refChord, ...
    '1.146-7');
AE_data{end + 1} = sprintf('$ MKAERO1: Mach = %.2f, list of 6 red frequencies', ...
    mach);
AE_data{end + 1} = sprintf('%-8s%-8.2f', 'MKAERO1', mach);
AE_data{end + 1} = sprintf('%-8s%-8.4f%-8.4f%-8.4f%-8.4f%-8.4f%-8.4f', ...
 ' ', redFreq(1), redFreq(2), redFreq(3), redFreq(4), redFreq(5), redFreq(6));
AE_data{end + 1} = '$ CAERO5: Piston Theory';
AE_data{end + 1} = sprintf('%-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s', '$', ...
    'EID', 'PID', 'CP', 'NSPAN', 'LSPAN', 'NTHRY', 'NTHICK');
AE_data{end + 1} = CAERO101{1};
AE_data{end + 1} = CAERO101{2};
AE_data{end + 1} = '';
AE_data{end + 1} = CAERO102{1};
AE_data{end + 1} = CAERO102{2};
AE_data{end + 1} = '';
AE_data{end + 1} = CAERO103{1};
AE_data{end + 1} = CAERO103{2};
AE_data{end + 1} = '';
AE_data{end + 1} = sprintf('%-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s', '$', ...
    'X1', 'Y1', 'Z1', 'X12', 'X4', 'Y4', 'Z4', 'X43');

AE_data{end + 1} = '';
AE_data{end + 1} = sprintf('%-8s%-8s%-8d', 'PARAM', 'LMODES', LMODES);
AE_data{end + 1} = sprintf('%-8s%-8d%-8s%-8s%-8d%-8s%-8s%-8s%-8s', 'EIGRL', 10, ...
    ' ', ' ', 18, ' ', ' ', ' ', 'MAX');

AE_data{end + 1} = sprintf('%-8s%-8s%-8s%-8s', '$', 'PID', 'NALPHA', 'LALPHA');
AE_data{end + 1} = sprintf('%-8s%-8d%-8d%-8d', 'PAERO5', 1, 1, 101);
AE_data{end + 1} = sprintf('%-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s', ' ', '0.', ...
    '0.', '0.','0.','0.' ,'0.','0.','0.');
AE_data{end + 1} = ['$       CAOC1   CAOC2   CAOC3   CAOC4   CAOC5   CAOC6   ' ...
    'CAOC7   CAOC8'];
AE_data{end + 1} = ' ';
AE_data{end + 1} = sprintf('%-8s%-8d%-8.2f%-8s', 'AEFACT', 101, mach, '0.');
AE_data = [AE_data, SET100TAB];
% thickness integrals:
AE_data{end + 1} = sprintf('%-8s%-8d%s', 'AEFACT', ...
    700, NTHICK1str);
AE_data{end + 1} = sprintf('%-8s%-8d%s', 'AEFACT', ...
    701, NTHICK2str);
AE_data{end + 1} = sprintf('%-8s%-8d%s', 'AEFACT', ...
    702, NTHICK3str);
AE_data{end + 1} = ' ';
AE_data{end + 1} = sprintf('%-8s%-8s%-8s%-8s%-8s%-8s%-8s', '$ ', 'ID', ...
    'CAEROID', 'BOX1', 'BOX2', 'SET ID');
AE_data{end + 1} = sprintf('%-8s%-8d%-8d%-8d%-8d%-8d%-8s', 'SPLINE1', ...
    101, 101, 101, 108, 100, "0.");
AE_data{end + 1} = sprintf('%-8s%-8d%-8d%-8d%-8d%-8d%-8s', 'SPLINE1', ...
    201, 201, 201, 208, 100, "0.");
AE_data{end + 1} = sprintf('%-8s%-8d%-8d%-8d%-8d%-8d%-8s', 'SPLINE1', ...
    301, 301, 301, 308, 100, "0.");
AE_data{end + 1} = ' ';
AE_data{end + 1} = sep;
AE_data{end + 1} = '$ FLUTTER SETUP';
AE_data{end + 1} = sep;
AE_data{end + 1} = sprintf('%-8s%-8s%-8s%-8s%-8s%-8s%-8s', '$', 'ID', 'METHOD', ...
    'DENSITY', 'MACH', 'VEL', 'INTERP');
AE_data{end + 1} = sprintf('%-8s%-8d%-8s%-8d%-8d%-8d%-8s', 'FLUTTER', 30, 'PK', ...
    1, 2, 3, 'L');
AE_data{end + 1} = ' ';
AE_data{end + 1} = sprintf('%-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s%-8s', '$', ...
    'SID', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7');
AE_data{end + 1} = sprintf('%-8s%-8d%-8.4f', 'FLFACT', 1, relDensity);
AE_data{end + 1} = sprintf('%-8s%-8d%-8.2f', 'FLFACT', 2, mach);

velTab = {};
velTab{end + 1} = sprintf('%-8s%-8d', 'FLFACT', 3);
k = 1;
for velIDX = 1:length(velocities)
    velTab{k} = [velTab{k}, sprintf('%-8.1f', ft2in * velocities(velIDX))];
    if mod(velIDX + 1, 8) == 0
        k = k + 1;
        velTab{k} = sprintf('%-8s', ' ');
    end
end
AE_data = [AE_data, velTab]';
AE_data{end + 1} = sprintf('%-8s%-8s%8.1f', 'PARAM', 'VREF', ft2in);
% OVERWRITING BDF DATA:
if formatted
    data = [data(1:comIDX - 1)
    AE_data
    oldDat(endIDX)];
else
    data = [data(1:end - 1)
    AE_data
    oldDat(endIDX)];
end
fid = fopen(bdfFile, 'w');
fprintf(fid, '%s\n', data{:});
fclose(fid);
end