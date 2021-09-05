% ---
% jupyter:
%   jupytext:
%     text_representation:
%       extension: .m
%       format_name: percent
%       format_version: '1.3'
%       jupytext_version: 1.4.0
%   kernelspec:
%     display_name: Matlab
%     language: matlab
%     name: matlab
% ---

% %% [markdown]
% # Figures from "Spatially organized multicellular immune hubs in human colorectal cancer"
% ## Code for key pannels included in Figure 2

% %%
code = '../code/';
codeExternal =  '../external';
data = '../data/';



% %%
run([code 'util/run_set_figure_defaults.m']);

% Other options (e.g. svg or eps) are possible:
outSuf = { '-dpng' } 
outDirPlot = { '../results/figures/' };

% %%
addpath(genpath(code))

% %%
addpath(genpath(codeExternal))

% %% [markdown]
% ## Loading main data object

% %%

colon10x_default = readDataRobj([ data 'colon10x_default/' ])


% %% [markdown]
% ## Loading tSNE coordinate files

% %%


[zFile,~,zFname] = fileList([ data 'cNMF_tSNE/*.gz']);

zFname = regexprep(zFname,'.tsv.gz','');
zFname = regexprep(zFname,'crc295v4_cNMF_tSNE_','');
zFname = regexprep(zFname,'allImm','Imm');
zFname{end} = 'global';

for zi = 1:length(zFname)
   [zRaw,~,zH] = fastMatRead(zFile{zi});    
   tSNE_coord.(zFname{zi}).sampleID = zH;
   tSNE_coord.(zFname{zi}).ydata = zRaw;
end

tSNE_coord

% %% [markdown]
% ## Load ccNMF summary file

% %%
ccNMFexpAlt = load([ data 'matlab/crc10x_c295v4_basic_ccNMFv6_reExp_qN_subSet.mat'])

zSubG = ismember(ccNMFexpAlt.ensgID,colon10x_default.ensgID);
ccNMFexpAlt = structSubSelectMat(ccNMFexpAlt,zSubG);
assert(isequal(ccNMFexpAlt.ensgID,colon10x_default.ensgID));


% %%

colon10x_default.metatable.MMRStatusTumor = colon10x_default.metatable.MMRStatus;
colon10x_default.metatable.MMRStatusTumor = regexprep(colon10x_default.metatable.MMRStatusTumor,'NA','Normal');

zBID = mergeStringPair(colon10x_default.metatable.MMRStatusTumor,colon10x_default.metatable.PatientTypeID);


% %%
%% Figre 2 
% Immune composition 
zOutP = 'Figure_2';

zSelImmB = { 'cB1' 'cB2' 'cB3'  'cP1' 'cP2' 'cP3' };
zSelImmM = { 'cMA01' 'cM01' 'cM02' 'cM03' 'cM04' 'cM05' 'cM06' 'cM07' 'cM08' 'cM09' 'cM10' };
zSelImmT = {'cTNI01' 'cTNI02' 'cTNI03' 'cTNI04' 'cTNI05' 'cTNI06' 'cTNI07' 'cTNI08' 'cTNI09' 'cTNI10' 'cTNI11' 'cTNI12' 'cTNI13' 'cTNI14' 'cTNI15' 'cTNI16' 'cTNI17' 'cTNI18' 'cTNI19' 'cTNI20' 'cTNI21' 'cTNI22' 'cTNI23' 'cTNI24' 'cTNI25' 'cTNI26'};

zSelImmAll = [ zSelImmT zSelImmB zSelImmM ];

%

zGroupT = groupTable(colon10x_default.metatable,'PatientTypeID');

%

zOrderCellAll = [ zGroupT.PatientTypeID(strcmp(zGroupT.MMRStatusTumor,'MMRd')); zGroupT.PatientTypeID(strcmp(zGroupT.MMRStatusTumor,'MMRp'));  zGroupT.PatientTypeID(strcmp(zGroupT.MMRStatusTumor,'Normal'))];

zOrderCell = [ zGroupT.PatientTypeID(strcmp(zGroupT.MMRStatusTumor,'MMRd')); zGroupT.PatientTypeID(strcmp(zGroupT.MMRStatusTumor,'MMRp'));];

% %% [markdown]
% ### Supplemental Figure 2B - Heatmap/dot comparing imune composition

% %%
%plot -s 1800,1800
zBsmp = colon10x_default.metatable.PatientTypeID;
[zB,~,~,zCnt] = fastUnique(zBsmp);
zIsOver1K = ismember(zBsmp,zB(zCnt>1000));

zSel = zIsOver1K & ~strgrep(cellstr(colon10x_default.annot.clTopLevel),'Epi|Strom',1);    
sum(zSel)

% tabFilter(colon10x_default.batchID,zSel)


% Comp vs normal 

zClust = regexprep(colon10x_default.annot.clFull(zSel),' .*','');
zPid = zBsmp(zSel);


zCC = crossTabTable(zClust,zPid,0,3);
% zSelN = strgrep(zCC.Properties.VariableNames,'N_'); !!
zSelN = strgrep(zCC.Properties.VariableNames,'_N');
zCC.Properties.RowNames;
zCCmat = table2array(zCC)+5; 

zCCnormal = zCCmat(:,zSelN);
zCCnormal = zCCnormal./sum(zCCnormal);

zCCnormal = median(zCCnormal,2);
zCCnormal = zCCnormal./sum(zCCnormal);

zCCexp = zCCnormal*sum(zCCmat);
zCCpr = (zCCmat - zCCexp)./sqrt(zCCexp);

zCCfreq = table2array(zCC);
zCCfreq = zCCfreq./sum(zCCfreq);


zLim = [ 0.05 0.1 0.15 0.2 0.25 ];

zCCcount = arrayfun(@(x)zCCfreq > x,zLim,'unif',0);


zCCfreqM = zCCfreq;
zCCfreqM = min(zCCfreqM,0.3);

zCCalt = array2table(zCCfreqM);
zCCalt.Properties.VariableNames = zCC.Properties.VariableNames;
zCCalt.Properties.RowNames = zCC.Properties.RowNames;

% 

[~,zia,zib] = intersect(zSelImmAll,zCC.Properties.RowNames);
zordT = zib(argsort(zia));

%

[~,zia,zib] = intersect(zOrderCell,zCC.Properties.VariableNames);
zCellOrder = zib(argsort(zia));

%%

% %%
%plot -s 1800,1800

zopts = [];
zopts.fontSize = 12;
zopts.pdistX = 'Euc';
zopts.pdistY = 'Spe';
zopts.pseudoCount = 0;
% zopts.colorType = 'pearson';
zopts.linkageX = 'ward';
zopts.linkageY = 'average';
zopts.showTxt = 0;

zopts.doSortX = 0;
zopts.doSortY = 0;
zopts.doLeafOptimalOrderX = 0;
zopts.doLeafOptimalOrderY = 0;

zopts.annotX.sampleID =  zGroupT.PatientTypeID;
zopts.annotX.MMR = zGroupT.MMRStatusTumor;
zopts.colorHard = [ -50 50 ];

zopts.extPR = zCCpr(zordT,zCellOrder);

zopts.colormap = redblue(128);
zopts.doDotPlot = 1;
zopts.dotPlotFactor = 500;
zopts.dotPlotColor = [0 0 0];
zopts.XTickLabelRotation = 90;


[zfig,zOutTable,zOrdX,zOrdY] = plot_crosstab_heatmap_alt([],[],[],zCCalt(zordT,zCellOrder),zopts);

%%

delete(zfig.Children(6:7));


% Legend

zxval = [ 0.01 0.05 0.1 0.15 0.2 0.25 0.3 ];
zxv = zxval*zopts.dotPlotFactor;
zn =length(zxv);
zxy = 1:zn;
zxx = ones(1,zn);

hold on;
zax = subplot('Position',[0.05 0.05 0.15 0.5]);
zaxp = scatter(zxx,zxy,zxv);
%
ylim([ 0 zn*3]);
xlim([ 0.75 1.25]);

%
zaxp.CData = [ 0 0 0];
zaxp.LineWidth = 1;

%
zvtxt = arrayfun(@(x)sprintf('%g',x),zxval,'unif',0);
zvtxt{end} = [ '>' zvtxt{end} ]

%
text(zxx+0.1,zxy,zvtxt,'fontsize',14);
axis off;

%% 

zOutPlot = sprintf('%s/SFig2b_ImmuneCellComp_enrichment_dot',zOutP);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);

% %% [markdown]
% ### Figure 2A Immune composition 

% %%
zClust = regexprep(colon10x_default.annot.clFull(zSel),' .*','');
zPid = zBsmp(zSel);

zCC = crossTabTable(zClust,zPid,0,3);
zSelN = strgrep(zCC.Properties.VariableNames,'_N');

zCCmat = table2array(zCC)+5; 

zCCnormal = zCCmat(:,zSelN);
zCCnormal = zCCnormal./sum(zCCnormal);

zCCnormal = median(zCCnormal,2);
zCCnormal = zCCnormal./sum(zCCnormal);

zCCexp = zCCnormal*sum(zCCmat);
zCCpr = (zCCmat - zCCexp)./sqrt(zCCexp);


%% 

[~,zia,zib] = intersect(zSelImmAll,zCC.Properties.RowNames);
zordT = zib(argsort(zia));


[~,zia,zib] = intersect(zOrderCell,zCC.Properties.VariableNames);
zCellOrder = zib(argsort(zia));

zCCfreq = table2array(zCC);
zCCfreq = zCCfreq./sum(zCCfreq);

%
zCCfreqM = zCCfreq;

zCCalt = array2table(zCCfreqM);
zCCalt.Properties.VariableNames = zCC.Properties.VariableNames;
zCCalt.Properties.RowNames = zCC.Properties.RowNames;

mmrMap = containers.Map(zGroupT.PatientTypeID,zGroupT.MMRStatusTumor);

%

% %%
zType = nanvalues(mmrMap,zCC.Properties.VariableNames);


zopts = [];
% zopts.aggrFunc = @(X,dim)quantile(X,0.5,dim);
zopts.aggrFunc = @(X,dim)mean(X,dim);
[zHfreq,~,zTypeN] = summarize_group(zCCfreq,zType,zopts);
max(zHfreq)

zopts = [];
zopts.aggrFunc = @(X,dim)mean(X,dim);
[zHpr,~,zTypeNpr] = summarize_group(zCCpr,zType,zopts);

% %%
zord = [ 3 2 1 ];

zTypeN = zTypeN(zord);

zHfreq = array2table(zHfreq(:,zord),'VariableNames',zTypeN);
% zHpr = array2table(zHpr(:,zord),'VariableNames',zTypeN);

zHpr = zHpr(:,zord);
zHpr(:,1) = 0; % Zero normal
zHfreq.Properties.RowNames = zCC.Properties.RowNames;

% zHpr.Properties.RowNames = zCC.Properties.RowNames;

[~,zia,zib] = intersect(zSelImmAll,zCC.Properties.RowNames);
zordT = zib(argsort(zia));

zopts = [];
zopts.fontSize = 12;
zopts.pdistX = 'Euc';
zopts.pdistY = 'Spe';
zopts.pseudoCount = 0;
% zopts.colorType = 'pearson';
zopts.linkageX = 'ward';
zopts.linkageY = 'average';
zopts.showTxt = 0;
zopts.doSortX = 0;
zopts.doSortY = 0;
zopts.doLeafOptimalOrderX = 0;
zopts.doLeafOptimalOrderY = 0;

zopts.colorHard = [ -50 50 ];

zopts.extPR = zHpr(zordT,:);

zopts.colormap = redblue(128);
zopts.doDotPlot = 1;
zopts.dotPlotFactor = 1000;
zopts.dotPlotColor = [0 0 0];
zopts.gridLines = '-';

[zfig,zOutTable,zOrdX,zOrdY] = plot_crosstab_heatmap_alt([],[],[],zHfreq(zordT,:),zopts);

%% 

delete(zfig.Children(3:5));
zfig.Position = [200 10 500 1700];

zfig.Children(2).Position(1) = 0.55;
zfig.Children(2).Position(2) = 0.1
zfig.Children(2).Position(3) = 0.23;
zfig.Children(2).Position(4) = 0.86

for zk = 2:length(zfig.Children(2).Children)-1
    zfig.Children(2).Children(zk).Color = [1 1 1];
end

% Legend
zxval = [ 0.01 0.05 0.1 0.15 0.2 0.25 0.3 ];
zxv = zxval*zopts.dotPlotFactor;
zn =length(zxv);
zxy = 1:zn;
zxx = ones(1,zn);

hold on;
zax = subplot('Position',[0.05 0.05 0.15 0.5]);
zaxp = scatter(zxx,zxy,zxv);
%
ylim([ 0 zn*3]);
xlim([ 0.75 1.25]);

%
zaxp.CData = [ 0 0 0];
zaxp.LineWidth = 1

%
zvtxt = arrayfun(@(x)sprintf('%g',x),zxval,'unif',0);
zvtxt{end} = [ '>' zvtxt{end} ];

%
text(zxx+0.1,zxy,zvtxt,'fontsize',14);
axis off;

%% 

zOutPlot = sprintf('%s/Fig2a_ImmuneCellComp_enrichment_mean_dot',zOutP);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);

% %% [markdown]
% ### Figure 2B -- Myeloid tSNE 

% %%
%plot -s 2820,620
zcf = 'Myeloid';
[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);

zYdata = tSNE_coord.(zcf).ydata(zib,:);
zCl = colon10x_default.annot.clFull(zia);
zMMRtype = colon10x_default.metatable.MMRStatusTumor(zia);

zfig = figure('Position',[0 0 2820 620]);
subplot1(1,4);

zopts = []
zopts.newPlot = 0;
zopts.pSize = 8;
zopts.doAlpha = 0.7;
zopts.doPosTxt = 1;
zopts.fontSize = 16;
zCmap = brewermap(luniq(zCl),'Dark2');

subplot1(1);
plot_tsne_scatter(zYdata,zCl,zCmap,[],zopts);

zyl = ylim();
zxl = xlim();

zopts.newPlot = 0;
zopts.pSize = 8;
zopts.doPosTxt = 0;
zopts.fontSize = 16;
legend off 

subplot1(2);
zSelS = strcmp(zMMRtype','Normal');
plot_tsne_scatter(zYdata(zSelS,:),zCl(zSelS),zCmap,[],zopts);
text(0.05,0.05,'Normal','sc');

xlim(zxl);
ylim(zyl);


subplot1(3);
zSelS = strcmp(zMMRtype,'MMRp');
plot_tsne_scatter(zYdata(zSelS,:),zCl(zSelS),zCmap,[],zopts);
text(0.05,0.05,'MMRp','sc');

legend off 
xlim(zxl);
ylim(zyl);

subplot1(4);
zSelS = strcmp(zMMRtype,'MMRd');
plot_tsne_scatter(zYdata(zSelS,:),zCl(zSelS),zCmap,[],zopts);
text(0.05,0.05,'MMRd','sc');

legend off 
xlim(zxl);
ylim(zyl);


%% 

zOutPlot = sprintf('%s/Fig2b_tSNE_%s_MMRsplit',zOutP,zcf);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);


% %% [markdown]
% ### Figure 2C -- Myeloid program compsition violin and tSNE

% %%
%plot -s 2600,1200
zcf = 'Myeloid'
[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);

zYdata = tSNE_coord.(zcf).ydata(zib,:);
zCl = colon10x_default.annot.clMidway(zia);
zMMRtype = colon10x_default.metatable.MMRStatusTumor(zia);
zMMRtype = colon10x_default.metatable.MMRStatusTumor(zia);

zHmat = ccNMFexpAlt.Hmat.(zcf)';
zHmat = zHmat(zia,:);

zWnames = ccNMFexpAlt.wNamesSt.(zcf);

% 

zBID = mergeStringPair(colon10x_default.metatable.MMRStatusTumor,colon10x_default.metatable.PatientTypeID);
zBIDsub = zBID(zia);

%

zopts = [];
zopts.aggrFunct = @(X,dim)quantile(X,0.75,dim);
[zHq75,~,zN] = summarize_group(zHmat',zBIDsub);
zHq75 = zHq75';

%

zNtype = regexprep(zN,'_.*','');

% 
% zopts = [];
% zopts.maxRow = 11;
% zopts.dodge = 0.6;
% plot_violin_simple(zWnames',zHq75,zNtype,zopts);
%

zfig = figure('Position',[ 0 0 2600 1200]);

zopts = [];
zopts.widthV = 0.9;
zopts.widthBox = 0.6;
zopts.doJitter = 0.4;
zopts.cmap = colorSet(zColSet.crcTypeC);
clear g 
for zi = 1:11 
    
    g(1,zi) = gramm('x',zNtype,'y',zHq75(:,zi),'Color',zNtype);
    g(1,zi).stat_violin2('width',zopts.widthV);
    g(1,zi).stat_boxplot('width',zopts.widthBox);
    g(1,zi).geom_jitter('width',zopts.doJitter);
    g(1,zi).set_names('x','','y',zWnames{zi});
end

for zi = 12:22
    zp = zi - 11;
    g(2,zp) = gramm('x',zNtype,'y',zHq75(:,zi),'Color',zNtype);
    g(2,zp).stat_violin2('width',zopts.widthV);
    g(2,zp).stat_boxplot('width',zopts.widthBox);
    g(2,zp).geom_jitter('width',zopts.doJitter);
    g(2,zp).set_names('x','','y',zWnames{zi});
end

g.set_color_options('map',zopts.cmap,'n_color',size(zopts.cmap,1));
g.set_layout_options('legend',false);

outG = g.draw();

%%

zfig = gcf();
zOutPlot = sprintf('%s/Fig2c_expProgramPidViolin_myeloid_MMRsplit',zOutP);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);
% close(zfig)

% %%
%plot -s 2400,2000

zcf = 'Myeloid';
zp = 'Myeloid';

[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);

zYdata = tSNE_coord.(zcf).ydata(zib,:);
zHval = ccNMFexpAlt.Hmat.(zp)';

zopts = [];
zopts.pSize = 7;
zopts.doAlpha = 0.7;
zopts.qTrimLimNNZ = 0.999;
zopts.titleText = ccNMFexpAlt.wNamesSt.(zp);
zopts.nCol = 4;
zopts.nRow = 3;
zopts.tiledLayout = 0;
zopts.plotSize = [ 1 1 2400 2000 ];
zopts.cmap = flipud(cptcmap('deep'));

zfig = plot_tsne_scatter_multi(zYdata,zHval(zia,:),[],[],zopts);

%% 

zOutPlot = sprintf('%s/Fig2c_tSNE_%s_Hmat_%s',zOutP,zcf,zp);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);
% cellfun(@(x)close(x),zfig)

% %% [markdown]
% ### Figure 2D -- T/NK/ILC tSNE 

% %%
%plot -s 2800,620
zcf = 'TNKILC'
[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);

zYdata = tSNE_coord.(zcf).ydata(zib,:);
zCl = colon10x_default.annot.clMidway(zia);
zMMRtype = colon10x_default.metatable.MMRStatusTumor(zia);

zfig = figure('Position',[0 0 2821 620]);
subplot1(1,4);

zopts = []
zopts.newPlot = 0;
zopts.pSize = 8;
zopts.doAlpha = 0.7;
zopts.doPosTxt = 1;
zopts.fontSize = 16;
zCmap = brewermap(luniq(zCl),'Dark2')

subplot1(1);
plot_tsne_scatter(zYdata,zCl,zCmap,[],zopts);

zyl = ylim();
zxl = xlim();

zopts.newPlot = 0;
zopts.pSize = 8;
zopts.doPosTxt = 0;
zopts.fontSize = 16;

subplot1(2);
zSelS = strcmp(zMMRtype','Normal');
plot_tsne_scatter(zYdata(zSelS,:),zCl(zSelS),zCmap,[],zopts);
text(0.05,0.05,'Normal','sc');

legend off 
xlim(zxl)
ylim(zyl);


subplot1(3);
zSelS = strcmp(zMMRtype,'MMRp');
plot_tsne_scatter(zYdata(zSelS,:),zCl(zSelS),zCmap,[],zopts);
text(0.05,0.05,'MMRp','sc');

legend off 
xlim(zxl)
ylim(zyl);

subplot1(4);
zSelS = strcmp(zMMRtype,'MMRd');
plot_tsne_scatter(zYdata(zSelS,:),zCl(zSelS),zCmap,[],zopts);
text(0.05,0.05,'MMRd','sc');

legend off 
xlim(zxl)
ylim(zyl);

%% 

zOutPlot = sprintf('%s/Fig2d_tSNE_%s_MMRsplit',zOutP,zcf);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);
% close(zfig)

% %%
%% Figure 2e - TODO 

% %% [markdown]
% ### Figure 2f

% %%
%plot -s 2600,1200
zcf = 'TNKILC'
zp = 'T';
[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);

zYdata = tSNE_coord.(zcf).ydata(zib,:);
zCl = colon10x_default.annot.clMidway(zia);
zMMRtype = colon10x_default.metatable.MMRStatusTumor(zia);
zMMRtype = colon10x_default.metatable.MMRStatusTumor(zia);

zHmat = ccNMFexpAlt.Hmat.(zp)';
zHmat = zHmat(zia,:);

zWnames = ccNMFexpAlt.wNamesSt.(zp);

%

zBID = mergeStringPair(colon10x_default.metatable.MMRStatusTumor,colon10x_default.metatable.PatientTypeID);
zBIDsub = zBID(zia);

%

zopts = [];
zopts.aggrFunct = @(X,dim)quantile(X,0.75,dim);
[zHq75,~,zN] = summarize_group(zHmat',zBIDsub);
zHq75 = zHq75';

%

zNtype = regexprep(zN,'_.*','');

%% 
% zopts = [];
% zopts.maxRow = 11;
% zopts.dodge = 0.6;
% plot_violin_simple(zWnames',zHq75,zNtype,zopts)

%

zfig = figure('Position',[ 0 0 2600 1200]);

zopts = [];
zopts.widthV = 0.9;
zopts.widthBox = 0.6;
zopts.doJitter = 0.4;
zopts.cmap = colorSet(zColSet.crcTypeC);
clear g 
for zi = 1:11 
    
    g(1,zi) = gramm('x',zNtype,'y',zHq75(:,zi),'Color',zNtype);
    g(1,zi).stat_violin2('width',zopts.widthV);
    g(1,zi).stat_boxplot('width',zopts.widthBox);
    g(1,zi).geom_jitter('width',zopts.doJitter);
    g(1,zi).set_names('x','','y',zWnames{zi});
end

for zi = 12:21
    zp = zi - 11;
    g(2,zp) = gramm('x',zNtype,'y',zHq75(:,zi),'Color',zNtype);
    g(2,zp).stat_violin2('width',zopts.widthV);
    g(2,zp).stat_boxplot('width',zopts.widthBox);
    g(2,zp).geom_jitter('width',zopts.doJitter);
    g(2,zp).set_names('x','','y',zWnames{zi});
end

g.set_color_options('map',zopts.cmap,'n_color',size(zopts.cmap,1));
g.set_layout_options('legend',false);

outG = g.draw();

%

zfig = gcf();
zOutPlot = sprintf('%s/Fig2f_expProgramPidViolin_TNKILC_MMRsplit',zOutP);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);
% close(zfig)


% %%
% Todo update colormap
%plot -s 2400,2000

zcf = 'TNKILC';
zp = 'T';

[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);

zYdata = tSNE_coord.(zcf).ydata(zib,:);
zHval = ccNMFexpAlt.Hmat.(zp)';

zopts = []
zopts.pSize = 7;
zopts.doAlpha = 0.7;
zopts.qTrimLimNNZ = 0.999;
zopts.titleText = ccNMFexpAlt.wNamesSt.(zp);
zopts.nCol = 4;
zopts.nRow = 3;
zopts.tiledLayout = 0;
zopts.plotSize = [ 1 1 2400 2000 ];
zopts.cmap = flipud(cptcmap('deep'));

zfig = plot_tsne_scatter_multi(zYdata,zHval(zia,:),[],[],zopts);

%

zOutPlot = sprintf('%s/Fig2f_tSNE_%s_Hmat_%s',zOutP,zcf,zp);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);
% cellfun(@(x)close(x),zfig)
