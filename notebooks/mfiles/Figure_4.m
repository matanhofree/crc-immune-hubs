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
% ## Code for key pannels included in Figure 4 -- Epithelial cells 

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
zOutP = 'Figure_4'
run run_set_figure_defaults

zMergeT = struct2table(groupTable(colon10x_default.metatable,'PatientTypeID'));

colon10x_default.metatable.MMRStatusTumor = colon10x_default.metatable.MMRStatus;
colon10x_default.metatable.MMRStatusTumor = regexprep(colon10x_default.metatable.MMRStatusTumor,'NA','Normal');

zBID = mergeStringPair(colon10x_default.metatable.MMRStatusTumor,colon10x_default.metatable.PatientTypeID);


% %% [markdown]
% ### Figure 4A

% %%
%plots -s 2400,2000

zcf = 'Epi'
[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);

zYdata = tSNE_coord.(zcf).ydata(zib,:);
zCl = colon10x_default.annot.clFull(zia);
zMMRtype = colon10x_default.metatable.MMRStatusTumor(zia);

zfig = figure('Position',[0 0 2400 2000])
subplot1(2,2);

zopts = []
zopts.newPlot = 0;
zopts.pSize = 8;
zopts.doAlpha = 0.7;
zopts.doPosTxt = 0;
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

zOutPlot = sprintf('%s/Fig4_tSNE_%s_MMRsplit',zOutP,zcf);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);
% close(zfig)

% %% [markdown]
% ### Figure 4A

% %%
%plot -s 2400,750
zcf = 'Epi'
[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);

zYdata = tSNE_coord.(zcf).ydata(zib,:);
zCl = colon10x_default.annot.clFull(zia);
zMMRtype = colon10x_default.metatable.MMRStatusTumor(zia);

zBIDsub = zBID(zia);

zClT = zCl;
zClT(~strcmp(zMMRtype,'Normal')) = {'Maligant'};

%%

zfig = figure('Position',[0 0 2400 750 ])

subplot1(1,3);

subplot1(1)

zopts = []
zopts.newPlot = 0;
zopts.pSize = 5;
zopts.doAlpha = 0.5;
zopts.doPosTxt = 0;
zopts.fontSize = 16;

zCmap = colorSet(zColSet.crcTypeC);

plot_tsne_scatter(zYdata,zMMRtype,zCmap,[],zopts);

subplot1(2)
zCmapRand = [ colorSet(zColSet.rainbow21); colorSet(zColSet.rainbow18); colorSet(zColSet.rainbow14);];
zCmapRand = [ zCmapRand; zCmapRand ];
zCmapRand = zCmapRand(randperm(length(zCmapRand)),:);

plot_tsne_scatter(zYdata,zBIDsub,zCmapRand,[],zopts);
legend off
text(0.05,0.05,'Color by patient sample','sc');

subplot1(3)
zCmap = brewermap(luniq(zClT),'Dark2')
plot_tsne_scatter(zYdata,zClT,zCmap,[],zopts);
text(0.05,0.05,'Color by cell-type','sc');

%%

zOutPlot = sprintf('%s/Fig4a_tSNE_%s_MMR_pid_cl',zOutP,zcf);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);
% close(zfig)

% %% [markdown]
% ### Figure 4B - Epithelial global program heatmap

% %%
%plot -s 2200,1800
[zB,~,~,zCnt] = fastUnique(colon10x_default.metatable.PatientTypeID);
zB(zCnt< 1000)
zSel1K = ismember(colon10x_default.metatable.PatientTypeID,setdiff(zB(zCnt > 1000),'C151_N'));


%%

zp = 'EpiTGlobalv5ForceK43';

zSel = zSel1K & strgrep(colon10x_default.annot.clTopLevel,'Epi');

zHmat = ccNMFexpAlt.Hmat.(zp)';
zHmat = zHmat(zSel,:);
zPidB = mergeStringPair(colon10x_default.metatable.MMRStatusTumor(zSel),colon10x_default.metatable.PatientTypeID(zSel));

zWnames = ccNMFexpAlt.wNamesSt.(zp);


zopts = [];
zopts.aggrFunc = @(X,dim)quantile(X,0.75,dim);
[zWq75mat,~,zWpid] = summarize_group(zHmat',zPidB,zopts);

%% 

zWnames = ccNMFexpAlt.wNamesSt.EpiTGlobalv5ForceK43;
zSelP = strgrep(zWpid,'_N');

zWpid(zSelP) = [];
zWq75mat(:,zSelP) = [];
%%

zopts = [] 
zopts.maxLabelY = 150;
zopts.maxLabelX = 150;
zopts.showTxt = 0;

zopts.annotX.sampleID = zMergeT.PatientTypeID;
zopts.annotX.MMR = zMergeT.MMRStatusTumor;

zopts.axisFontSize = 9;
zopts.colorLim = [ 0.02 0.98 ];

zWpidT = regexprep(zWpid,'^[^_]*_','');
%%

zX = nanzscore(zWq75mat,[],2);

%%
% TODO add adjacent normal

[zfigOut,outMat,outMat.zOrdTX,outMat.zOrdTY] = plot_heatmap_annot(zX,zWnames,zWpidT,[],[],zopts);

%%

zOutPlot = sprintf('%s/Fig4b_pEpiGlobal_q75_HC',zOutP);    
cFname = print_plot(zfigOut,zOutPlot,outDirPlot,outSuf{1},1);
% close(zfigOut)

% %% [markdown]
% ### Figure 4C -- Composition of inferred epithelial cell types

% %%
%plot -s 1800,600

zClEpi = colon10x_default.annot.clFull(zSel);
zPidEpi =colon10x_default.metatable.PatientTypeID(zSel);

% plot_crosstab_bar(zPidEpi,zClEpi,[],zopts)

zClPidTable = crossTabTable(zClEpi,zPidEpi,1,3);

%% 

[~,zia,zib] = intersect(zClPidTable.Properties.VariableNames,outMat.xNameOrd);

zord = zia(argsort(zib));

zClPidTorder = zClPidTable(:,zord);

zSelN = strgrep(zClPidTable.Properties.VariableNames,'_N');
zTabN = zClPidTable(:,zSelN);
zFreqN = median(table2array(zTabN),2);
zFreqN = zFreqN/sum(zFreqN);

zClPidTorder.Adj_Normal = zFreqN;

%% 

zfig = figure('Position',[ 1 1 1800 600]);
zax = bar(table2array(zClPidTorder)','stacked');
xName = zClPidTorder.Properties.VariableNames;
box off 
axis tight 
%%

zax = zfig.Children;
set(zax,'xtick',1:length(xName));
set(zax,'xticklabels', regexprep(xName,'_','-'));

zax.XTickLabelRotation = 60
legend(unique(zClEpi),'Location','northeastoutside','box','off','interpreter','none')

%%

zOutPlot = sprintf('%s/Fig4c_maligEpi_cellcomp_bar',zOutP);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);
% close(zfigOut)

% %% [markdown]
% ### Epithelial global program violin plots and tSNE

% %%
%plot -s 2600,1600

zcf = 'Epi';
zp = 'EpiTGlobalv5ForceK43';

[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);

zYdata = tSNE_coord.(zcf).ydata(zib,:);
zCl = colon10x_default.annot.clMidway(zia);
zMMRtype = colon10x_default.metatable.MMRStatusTumor(zia);

zHmat = ccNMFexpAlt.Hmat.(zp)';
zHmat = zHmat(zia,:);

zWnames = ccNMFexpAlt.wNamesSt.(zp);

%% 

zBID = mergeStringPair(colon10x_default.metatable.MMRStatusTumor,colon10x_default.metatable.PatientTypeID);
zBIDsub = zBID(zia);

%%

zopts = [];
zopts.aggrFunc = @(X,dim)quantile(X,0.75,dim);
[zHq75,~,zN] = summarize_group(zHmat',zBIDsub);
zHq75 = zHq75';

%%

zNtype = regexprep(zN,'_.*','');

%% 

zPn = luniq(zWnames);
zNrow = 4;
zNcol = ceil(zPn/zNrow);

%%


zfig = figure('Position',[ 0 0 2600 500*zNrow]);

zopts = [];
zopts.widthV = 0.9;
zopts.widthBox = 0.6;
zopts.doJitter = 0.4;
zopts.cmap = colorSet(zColSet.crcTypeC);
clear g 

for zi = 1:zPn
    zx = floor((zi-1)/zNcol) + 1;
    zy = mod(zi-1,zNcol) + 1;
    
    g(zx,zy) = gramm('x',zNtype,'y',zHq75(:,zi),'Color',zNtype);
    g(zx,zy).stat_violin2('width',zopts.widthV);
    g(zx,zy).stat_boxplot('width',zopts.widthBox);
    g(zx,zy).geom_jitter('width',zopts.doJitter);
    g(zx,zy).set_names('x','','y',zWnames{zi});
end


g.set_color_options('map',zopts.cmap,'n_color',size(zopts.cmap,1));
g.set_layout_options('legend',false);

outG = g.draw();

%%

zfig = gcf();
zOutPlot = sprintf('%s/Fig4d_expProgramPidViolin_Epi_MMRsplit',zOutP);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);
% close(zfig)

% %%
%plot -s 1900,1150

zcf = 'Epi';
zp = 'EpiTGlobalv5ForceK43';

[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);

zYdata = tSNE_coord.(zcf).ydata(zib,:);
zHval = ccNMFexpAlt.Hmat.(zp)';

zopts = []
zopts.pSize = 7;
zopts.doAlpha = 0.7;
zopts.qTrimLimNNZ = 0.999;
zopts.titleText = ccNMFexpAlt.wNamesSt.(zp);
zopts.nCol = 5;
zopts.nRow = 3;
zopts.tiledLayout = 0;
zopts.plotSize = [ 1 1 1900 1150 ];
zopts.cmap = flipud(cptcmap('deep'));

zfig = plot_tsne_scatter_multi(zYdata,zHval(zia,:),[],[],zopts);

%% 

zOutPlot = sprintf('%s/Fig4d_tSNE_%s_Hmat_%s',zOutP,zcf,zp);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);

% %%
