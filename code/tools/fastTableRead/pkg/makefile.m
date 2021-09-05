%% Make for fastTable Read function 

if ismac()

    mex -v fastTableReadV2.cpp -largeArrayDims -I/usr/local/include/ -lboost_iostreams
elseif isunix()
    
    mex -v -lz -largeArrayDims -I/ahg/regevdata/projects/sawyersProstateSC/extra/library/boost_1_76_0 /ahg/regevdata/projects/sawyersProstateSC/extra/library/boost_1_76_0/stage/lib/libboost_iostreams.a fastTableReadV3.cpp
end

%%
fname = '/ahg/regevdata/users/mhofree/local/toolbox-matlab/fastTableRead/test_data.small.txt';
[zz,zerr] = fastTxtRead(fname,'\t');

%% 

fname = './test/test_data.small.txt';
[zdata,zerr] = fastTableReadV3(2,fname,[],char(9));

%%

fname = './test/test_data.small.txt';
[zdata,zerr,zRowH,zColH] = fastTableReadV3(2,fname,[],char(9),0,1,1)

%%

fname = './test/test_data_sm.txt.gz';
[zOutD.zdata,zerr,zOutD.zRowH,zOutD.ColH] = fastTableReadV3(2,fname,[],char(9),0,1,1)

%%

fname = './test/test_data_sm.txt.gz';
[zOutSp.zdata,zerr,zOutSp.zRowH,zOutSp.zColH] = fastTableReadV3(3,fname,[],char(9),0,1,1);
zOutSp.zdata = zOutSp.zdata'

%% 

[zi,zj,zv] = find(zOutSp.zdata)

%%

[zx.zi,zx.zj,zx.zv] = find(zOutD.zdata');
%
struct2table(zx)

%% 

isequal(zOutD.zdata,full(zOutSp.zdata))
isequal(zOutD.zRowH,full(zOutSp.zRowH))
isequal(zOutD.ColH,full(zOutSp.zColH))

