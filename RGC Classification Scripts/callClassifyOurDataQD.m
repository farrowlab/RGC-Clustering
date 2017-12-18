%% Initiate
clear
clc

addpath('/media/areca_raid/LabCode/Matlab/Clustering')
addpath('/media/areca_raid/LabCode/Matlab/Clustering/fastAPv2')
addpath('/media/areca_raid/LabCode/Matlab/Clustering/FUZZCLUST')
addpath('/media/areca_raid/LabCode/Matlab/Clustering/minCEntropy')
addpath('/media/areca_raid/LabCode/Matlab/Clustering/spasm')

%% Adjustable variables
% pathes
dname = '/media/areca_raid/LabPapers/SCRouter/Data/dendriticArbors';
dbase = '/media/areca_raid/LabPapers/SCRouter/Data/database';
path2save = '/media/areca_raid/LabPapers/SCRouter/Data/clusteringResults';

date = '20171128';
%----- User Needs to Set These Variable -----%0=optu
types2analyse = [1 3]; % 1=PBg, 2=LP, 3=floxLP
includeMaybes = 1;
invertBistrat = 1;
conditions = 20;
zProfileThreshold = 0;
xSize = 28;
downZ = {0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0};
nz2 = 12; %12 24 30 40 60 120
% doDimReduction = {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0};
doDimReduction = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
dimreductionMeth = 'PCA'; % {'PCA', 'SVD', 'sPCA', 'tSNE', 'NNMF'}
NumberDimensions = 5;
CMeths = {'E-Linkage','E-Linkage','E-Linkage','E-Linkage',...
    'EM-GMM','EM-GMM','EM-GMM','EM-GMM',...
    'E-Linkage','E-Linkage','E-Linkage','E-Linkage',...
    'E-Linkage','E-Linkage','E-Linkage','E-Linkage',...
    'EM-GMM','EM-GMM','EM-GMM','EM-GMM','E-Linkage'};%, 'EM-GMM'};%, 'K-Medoids'}; % 'Cluster Methods: E-Linkage, K-Means, K-Medoids, Affinity-Propogation, Fuzzy -C, Fuzzy - GK
dist = { 'euclidean', 'euclidean', 'euclidean', 'euclidean',...
    '','','','', 'euclidean', 'euclidean', 'euclidean', 'euclidean',...
    'mahalanobis','mahalanobis','mahalanobis','mahalanobis',...
    '','','','','mahalanobis'}; % 'sqeuclidean', 'cityblock', 'correlation', 'cosine'}; {'mahalanobis', 'euclidean', 'sqEuclidean', 'correlation'

arbourType = 'arborDensity2';
arbourSpec = {'','','_tree','_tree','','','_tree','_tree',...
    '','','_tree','_tree','','','_tree','_tree','','','_tree','_tree',''};
clc
for c = 13%[1:4 9:20]
     if iscell(invertBistrat)
        azProfileThreshold = zProfileThreshold {c};
    else
        azProfileThreshold = zProfileThreshold;
    end
    if iscell(invertBistrat)
        ainvertBistrat = invertBistrat {c};
    else
        ainvertBistrat = invertBistrat;
    end
     if iscell(types2analyse)
        atypes2analyse = types2analyse {c};
    else
        atypes2analyse = types2analyse;
    end
    if iscell(includeMaybes)
        aincludeMaybes = includeMaybes {c};
    else
        aincludeMaybes = includeMaybes;
    end
     if iscell(downZ)
        adownZ = downZ {c};
    else
        adownZ = downZ;
    end
    if iscell(nz2)
        anz2 = nz2 {c};
    else
        anz2 = nz2; 
    end
    if iscell(doDimReduction)
        adoDimReduction = doDimReduction {c};
    else
        adoDimReduction = doDimReduction;
    end
    if iscell(dimreductionMeth)
        adimreductionMeth = dimreductionMeth {c};
    else
        adimreductionMeth = dimreductionMeth;
    end
    if iscell(NumberDimensions)
        aNumberDimensions = NumberDimensions {c};
    else
        aNumberDimensions = NumberDimensions;
    end
    if iscell(CMeths)
        clear aCMeths
        aCMeths{1,1} = CMeths {c};
    else
        aCMeths{1,1} = CMeths;
    end
    if iscell(dist)
        adist = dist {c};
    else
        adist = dist;
    end
    if iscell(arbourType)
        aarbourType = arbourType {c};
    else
        aarbourType = arbourType;
    end
    if iscell(xSize)
        axSize = xSize {c};
    else
        axSize = xSize;
    end
    if iscell(arbourSpec)
        aarbourSpec = arbourSpec {c};
    else
        aarbourSpec = arbourSpec;
    end
    arbourName = [aarbourType,aarbourSpec];
    
    if adownZ == 1
        part1 = ['d',int2str(anz2)];
    else
        part1 = ['d0'];
    end
    if adoDimReduction == 1
        part2 = ['r',adimreductionMeth,int2str(aNumberDimensions)];
    else
        part2 = ['r0'];
    end
    part3 = [aCMeths{1,1},adist];
    tmp = num2str(azProfileThreshold);
    if length(tmp)>1
        tmp(2) = 'p';
    end
    part4 = ['t',tmp];
    part5 = ['m',int2str(aincludeMaybes)];
    
    
    name2save = [date,'_',part1,'_',part2,'_',part3,'_',arbourName,'_',part4,'_',part5];
    disp([int2str(c),' | ',name2save])
    
    classifyRGCarborDensitiesQD(axSize,ainvertBistrat, aincludeMaybes,atypes2analyse,dname, dbase,path2save,adownZ,anz2,adoDimReduction,...
        adimreductionMeth,aNumberDimensions,aCMeths,adist,arbourName,name2save,azProfileThreshold)
    
end
