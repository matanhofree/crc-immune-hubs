%% Colorset

zColSet.rainbow14={'#882E72', '#B178A6', '#D6C1DE', '#1965B0', '#5289C7', '#7BAFDE', '#4EB265', '#90C987', '#CAE0AB', '#F7EE55', '#F6C141', '#F1932D', '#E8601C', '#DC050C'};
zColSet.rainbow15={'#114477', '#4477AA', '#77AADD', '#117755', '#44AA88', '#99CCBB', '#777711', '#AAAA44', '#DDDD77', '#771111', '#AA4444', '#DD7777', '#771144', '#AA4477', '#DD77AA'};
zColSet.rainbow18={'#771155', '#AA4488', '#CC99BB', '#114477', '#4477AA', '#77AADD', '#117777', '#44AAAA', '#77CCCC', '#777711', '#AAAA44', '#DDDD77', '#774411', '#AA7744', '#DDAA77', '#771122', '#AA4455', '#DD7788'};
zColSet.rainbow21={'#771155', '#AA4488', '#CC99BB', '#114477', '#4477AA', '#77AADD', '#117777', '#44AAAA', '#77CCCC', '#117744', '#44AA77', '#88CCAA', '#777711', '#AAAA44', '#DDDD77', '#774411', '#AA7744', '#DDAA77', '#771122', '#AA4455', '#DD7788'};
zColSet.crcTypeColor = {'#0072B2' '#B2DF8A' '#E41A1C' '#FB9A99' };
zColSet.crcTypeC = { '#E41A1C' '#0072B2' '#B2DF8A' };

zColSet.crcTypeColorName = {'MSS' 'Normal' 'MSI_MLH1Meth' 'MSI_MLH1NoMeth' };
zColSet.crcTypeColorName = {'MMRP' 'Normal' 'MMRd' 'MMRd_MLH1NoMeth' };


colorSet = @(zCol)cell2mat(cellfun(@(x)hexColor(x(2:end)),zCol,'uniformoutput',0))';
fixNames = @(x)regexprep(x,'_','-');

%%

zFontT = 16;
set(0, 'DefaultFigureColor', 'White', ...
'DefaultFigurePaperType', 'a4letter', ...
'DefaultAxesColor', 'white', ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', zFontT, ...
'DefaultAxesFontName', 'Ariel', ...
'DefaultAxesGridLineStyle', ':', ...
'DefaultAxesInterruptible', 'on', ...
'DefaultAxesLayer', 'Bottom', ...
'DefaultAxesNextPlot', 'replace', ...
'DefaultAxesUnits', 'normalized', ...
'DefaultAxesXcolor', [0, 0, 0], ...
'DefaultAxesYcolor', [0, 0, 0], ...
'DefaultAxesZcolor', [0, 0, 0], ...
'DefaultAxesVisible', 'on', ...
'DefaultLineColor', 'Red', ...
'DefaultLineLineStyle', '-', ...
'DefaultLineLineWidth', 1, ...
'DefaultLineMarker', 'none', ...
'DefaultLineMarkerSize', 8, ...
'DefaultTextColor', [0, 0, 0], ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize',  zFontT, ...
'DefaultTextFontName', 'Ariel', ...
'DefaultTextVerticalAlignment', 'middle', ...
'DefaultTextHorizontalAlignment', 'left');

%%

% cmap = cptcmap('deep')