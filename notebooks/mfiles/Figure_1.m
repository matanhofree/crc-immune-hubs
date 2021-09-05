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
% ## Code for key pannels included in Figure 1
%

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
%% Figure 1 

zOutP = 'Figure_1'

zMergeT = struct2table(groupTable(colon10x_default.metatable,'PatientTypeID'));
[zB,~,~,zCnt,zPos] = fastUnique(zMergeT.MMRStatus);

zBlues = brewermap(zCnt(1),'blues');
zGreens = brewermap(zCnt(2),'YlGn');
zReds = brewermap(zCnt(3),'OrRd');

zCombMap = [ zReds; zBlues; zGreens; ];

colon10x_default.metatable.MMRStatusTumor = colon10x_default.metatable.MMRStatus;
colon10x_default.metatable.MMRStatusTumor = regexprep(colon10x_default.metatable.MMRStatusTumor,'NA','Normal');
zBID = mergeStringPair(colon10x_default.metatable.MMRStatusTumor,colon10x_default.metatable.PatientTypeID);


% %%
%plot -s 1200,1000

zcf = 'global';
[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);

zYdata = tSNE_coord.(zcf).ydata(zib,:);

zCmap = brewermap(7,'Set1');

zopts = [];
zopts.pSize = 5;
zopts.doAlpha = 0.7;

plot_tsne_scatter(zYdata,colon10x_default.annot.clTopLevel(zia),zCmap,[],zopts);
zfig = gcf;

axis off 
axis tight

% Todo: fix colors

zOutPlot = sprintf('%s/Fig1a_tSNE_%s_topLevelCl',zOutP,zcf);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf,1);

% %%
zcf = 'global'
[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);

zYdata = tSNE_coord.(zcf).ydata(zib,:);

zCmap = colorSet(zColSet.crcTypeC);

zopts = [];
zopts.pSize = 5;
zopts.doAlpha = 0.7;

plot_tsne_scatter(zYdata,colon10x_default.metatable.MMRStatusTumor(zia),zCmap,[],zopts);
zfig = gcf;

axis off 
axis tight
% legend off

%% 
zOutPlot = sprintf('%s/Fig1b_tSNE_%s_MMR',zOutP,zcf);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);
% close(zfig)

% %%
%% MMR plot 

zcf = 'global'
[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);

zYdata = tSNE_coord.(zcf).ydata(zib,:);

zCmapRand = [ colorSet(zColSet.rainbow21); colorSet(zColSet.rainbow18); colorSet(zColSet.rainbow14);];
zCmapRand = zCmapRand(randperm(length(zCmapRand)),:);

zopts = [];
zopts.pSize = 5;
zopts.doAlpha = 0.7;

plot_tsne_scatter(zYdata,zBID(zia),zCmap,[],zopts);
zfig = gcf;

axis off 
axis tight
legend off

zOutPlot = sprintf('%s/Fig1c_tSNE_%s_byPidColor',zOutP,zcf);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);
% close(zfig)



% %%
%% Global by PID 
% MMRd, MMRp, and Normal patients in shaded ins reds, blues, and greens respectively. 


zcf = 'global'
[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);


zYdata = tSNE_coord.(zcf).ydata(zib,:);

zopts = [];
zopts.pSize = 7;
zopts.doAlpha = 0.7;

plot_tsne_scatter(zYdata,zBID(zia),zCombMap,[],zopts);
zfig = gcf;

axis off 
axis tight
legend off

%%

zOutPlot = sprintf('%s/tSNE_%s_byPid',zOutP,zcf);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);
% close(zfig)

% %%
%plot -s 2000,1500
%% Plot program NMF on tSNE
% Todo update colormap

zcf = 'global';
zp = 'EpiTMSIv4ForceK29';

[zia,zib] = comember(colon10x_default.sampleID,tSNE_coord.(zcf).sampleID);

zYdata = tSNE_coord.(zcf).ydata(zib,:);
zHval = ccNMFexpAlt.Hmat.(zp)';

zopts = [];
zopts.pSize = 7;
zopts.doAlpha = 0.7;
zopts.qTrimLimNNZ = 0.999;
zopts.titleText = ccNMFexpAlt.wNamesSt.(zp);

zopts.cmap = flipud(cptcmap('deep'));

zfig = plot_tsne_scatter_multi(zYdata,zHval(zia,:),[],[],zopts);

%% 

zOutPlot = sprintf('%s/Fig1c_tSNE_%s_Hmat_%s',zOutP,zcf,zp);    
cFname = print_plot(zfig,zOutPlot,outDirPlot,outSuf{1},1);
% cellfun(@(x)close(x),zfig)
