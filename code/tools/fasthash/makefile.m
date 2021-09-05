%%

mex -v fasthashStrPtr.c ... 
    -I./uthash ... 
    -largeArrayDims ...
    CFLAGS="-std=c99 -O2 -fPIC -fno-omit-frame-pointer -m64"

%% 
mex -v fasthashStr.c ... 
    -I./uthash ... 
    -largeArrayDims ...
    CFLAGS="-std=c99 -O2 -fPIC -fno-omit-frame-pointer -m64"

%%
mex -v fasthashInt.c ... 
    -I./uthash ... 
    -largeArrayDims ...
    CFLAGS="-std=c99 -O2 -fPIC -fno-omit-frame-pointer -m64"

%%
mex -v fastUniqueInt.cpp ...  
    -I/ahg/regevdata/users/mhofree/local/external-matlab/mexplus/include ...
    -largeArrayDims 

%%

mex -v fastUniqueIntFlat.cpp ...  
    -I/ahg/regevdata/users/mhofree/local/external-matlab/mexplus/include ...
    -largeArrayDims 
    
%%

%%

zz = randi(10000,100000,1);
zz = uint64(zz);

%%
tic;
[b, m, n, cnt,dup_list] = fastUniqueInt(zz);
toc
%%
[zb,zm,zn,zcnt,zdup_list] = uniquec(zz);


%%

isequal(zb,double(b'))
isequal(zcnt,double(cnt'))    
all(cellfun(@(x,y)isequal(x,double(y)),zdup_list,dup_list'))
%%
isequal(double(m),zm)

%%

%%
tic;
[b, m, n, cnt,dupListIdx,dupVect] = fastUniqueIntFlat(zz);
toc
%%

isequal(zb,double(b'))
isequal(zcnt,double(cnt'))    
isequal(zm,m)
isequal(zn,n)

%%

dup_list_re = arrayfun(@(zi)dupVect(dupListIdx(zi):dupListIdx(zi+1)-1)',1:length(dupListIdx)-1,'uniformoutput',0);
dup_list_re{end+1} = dupVect(dupListIdx(end):end)';
%%
all(cellfun(@(x,y)isequal(x,double(y)),zdup_list,dup_list_re))

%%


mex -v fastUniqueIntLL.cpp ...  
    -I/ahg/regevdata/users/mhofree/local/external-matlab/mexplus/include ...
    -largeArrayDims 


%%

tic;
[b, m, n, cnt,dupVect] = fastUniqueIntLL(zz);
toc
%%

isequal(zb,double(b))
isequal(zcnt,double(cnt))    
isequal(zm,m)
isequal(zn,n)

%%
dupListIdx = [ 1; cumsum(cnt)+1 ];
dup_list_re = arrayfun(@(zi)dupVect(dupListIdx(zi):dupListIdx(zi+1)-1)',1:length(dupListIdx)-1,'uniformoutput',0);

%%
% dup_list_re = arrayfun(@(zi)dupVect(dupListIdx(zi):dupListIdx(zi+1)-1)',1:length(dupListIdx)-1,'uniformoutput',0);
% dup_list_re{end+1} = dupVect(dupListIdx(end):end)';
%%
all(cellfun(@(x,y)isequal(x,double(y)),zdup_list,dup_list_re))


%%


mex -v fastUniqueIntMultiPass.cpp ...  
    -I/ahg/regevdata/users/mhofree/local/external-matlab/mexplus/include ...
    -largeArrayDims 


%%

tic;
[b, m, n, cnt,dupVect] = fastUniqueIntMultiPass(zz);
toc
%%

isequal(zb,double(b'))
isequal(zcnt,double(cnt'))    
isequal(zm,m')
isequal(zn,n')

%%
dupListIdx = [ 1; cumsum(cnt)+1 ];
dup_list_re = arrayfun(@(zi)dupVect(dupListIdx(zi):dupListIdx(zi+1)-1)',1:length(dupListIdx)-1,'uniformoutput',0);
%
all(cellfun(@(x,y)isequal(x,double(y)),zdup_list,dup_list_re))



%%
mex -v fastUniqueStr.cpp ... 
    -I/ahg/regevdata/users/mhofree/local/toolbox-matlab/fasthash/uthash ... 
    -I/ahg/regevdata/users/mhofree/local/external-matlab/mexplus/include ...
    -largeArrayDims 
    % CXXFLAGS="-fno-common -arch x86_64 -mmacosx-version-min=10.11 -fexceptions -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -fobjc-arc -std=c++11 -stdlib=libc++ -O2 -fPIC -fno-omit-frame-pointer -m64"
    

%%

load ('/ahg/regevdata/users/mhofree/projects/cancer_SC/results/2016_06_23_colon_new/temp/zGumi_list.mat')
%%

[b, m, n, cnt,dup_list] = fastUnique(zGumi);

%%


mex -v fasthashIntSimple.cpp ...  
    -I/ahg/regevdata/users/mhofree/local/external-matlab/mexplus/include ...
    -largeArrayDims 
    

%%

zKey = 1:100;
zz = randi(110,100000,1);
zz = uint64(zz);

zPos = fasthashIntSimple(zKey,zz);

%%

mex -v fasthashInt.c ... 
    -I/ahg/regevdata/users/mhofree/local/toolbox-matlab/fasthash/uthash ... 
    -largeArrayDims % ...
%    CFLAGS="-std=c99 -O2 -fPIC -fno-omit-frame-pointer -m64"

%%


zKey = 1:100;
zz = randi(100,100000,1);
zz = uint64(zz);

zPos = fasthashInt(uint64(zKey),zz);


%%

mex -v fastUniqueIntHash.c ... 
    -I/ahg/regevdata/users/mhofree/local/toolbox-matlab/fasthash/uthash ... 
    -largeArrayDims  % ...
   % CFLAGS="-std=c99 -O2 -fPIC -fno-omit-frame-pointer -m64 -fno-common -arch x86_64 -mmacosx-version-min=10.11 -fexceptions -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk"


%%
zz = randi(100,1000,1);
zz = uint64(zz);

%%
tic;
[b, m, n, cnt,dup_list] = fastUniqueInt(zz);
toc
%%
[zb,zm,zn,zcnt,zdup_list] = uniquec(zz);

%%



tic;
[b, m, n, cnt,dupVect] = fastUniqueIntHash(zz,length(unique(zz)));
toc
%%

isequal(zb,double(b'))
isequal(zcnt,double(cnt'))    
isequal(zm,m')
isequal(zn,n')

%%
dupListIdx = [ 1; cumsum(cnt)+1 ];
dup_list_re = arrayfun(@(zi)dupVect(dupListIdx(zi):dupListIdx(zi+1)-1)',1:length(dupListIdx)-1,'uniformoutput',0);
%
all(cellfun(@(x,y)isequal(x,double(y)),zdup_list,dup_list_re))

%% 

zz = randi(10000,1000000,1);
zz = uint64(zz);

%%
tic;
[b, m, n, cnt,dup_list,dupVect] = fastUnique(zz);
toc
%% Mac

mex -v fasthashInt.c ... 
    -I/Users/mhofree/local/external-matlab-sync/uthash/include ... 
    -largeArrayDims ...
    CFLAGS="-std=c99 -O2 -fPIC -fno-omit-frame-pointer -m64"

%%

mex -v fasthashStrPtr.c ... 
    -Iuthash/include ... 
    -largeArrayDims ...
    CFLAGS="-std=c99 -O2 -fPIC -fno-omit-frame-pointer -m64"

%%

mex -v fasthashInt.c -I/ahg/regevdata/users/mhofree/local/external-matlab/uthash/include ... 
    -largeArrayDims ...
    CFLAGS="-std=c99 -O2 -fPIC -fno-omit-frame-pointer -m64"

mex -v fasthashStrPtr.c -I/ahg/regevdata/users/mhofree/local/external-matlab/uthash/include ... 
    -largeArrayDims ...
    CFLAGS="-std=c99 -O2 -fPIC -fno-omit-frame-pointer -m64"

%% FastUniq OS catalina 

mex -v fastUniqueStr.cpp ...  
    -I/Users/mhofree/local/external-matlab-sync/mexplus/include ...
    -largeArrayDims 
    
%%
mex -v fastUniqueIntHash.c ...  
    -I/Users/mhofree/local/external-matlab-sync/mexplus/include ...
    -I/Users/mhofree/local/external-matlab-sync/uthash/include ...
    -largeArrayDims 
    