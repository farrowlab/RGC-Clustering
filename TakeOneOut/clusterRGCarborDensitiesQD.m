function     clusterRGCarborDensitiesQD(xSize,invertBistrat, includeMaybes,types2analyse,dname, dbase,path2save,downZ,nz2,doDimReduction,dimreductionMeth,NumberDimensions,CMeths,dist,arbourName,name2save,zProfileThreshold)
close all  
K = 21;
%% Load variables and find needed data
umPerPix = 21;
T1 = zProfileThreshold;

load(fullfile(dbase, 'usable.mat'));
usable = usable{1,2};
if includeMaybes == 1
    totake = find(usable>0);
else
    totake = find(usable==1);%absolute ids
end
load(fullfile(dbase, 'experiment.mat'));
load(fullfile(dbase, 'mouse.mat'));
load(fullfile(dbase, 'retina.mat'));
load(fullfile(dbase, 'scan.mat'));
load(fullfile(dbase, 'cellNumber.mat'));
load(fullfile(dbase,'areaHull.mat'));
areas = areaHull{1,2}; areas = areas(totake);
load(fullfile(dbase, 'type.mat'));
types = type{1,2}; types = types(totake);
load(fullfile(dbase,'zpeak1.mat'));
zpeaks1 = zpeak1{1,2}; zpeaks1 = zpeaks1(totake);
load(fullfile(dbase,'zpeak2.mat'));
zpeaks2 = zpeak2{1,2}; zpeaks2 = zpeaks2(totake);
load(fullfile(dbase,'bistratType.mat'));
bistratType = bistratType{1,2}; bistratType = bistratType(totake);
load(fullfile(dbase,'density.mat'));
density = density{1,2}; density = density(totake);
load(fullfile(dbase,'DSI.mat'));
DSI = DSI{1,2}; DSI = DSI(totake);
load(fullfile(dbase,'OSI.mat'));
OSI = OSI{1,2}; OSI = OSI(totake);
load(fullfile(dbase,'xCoordinates.mat'));
xCoordinates = xCoordinates{1,2}; xCoordinates = xCoordinates(totake);
load(fullfile(dbase,'yCoordinates.mat'));
yCoordinates = yCoordinates{1,2}; yCoordinates = yCoordinates(totake);

%% Get data

    %---------- cells to analyse ----------%
    cells2use = []; cells2useAbs = [];
    for t = types2analyse
        cells2use = [cells2use; find(types == t)];
        cells2useAbs = [cells2useAbs; totake(types == t)];
    end
    disp(['Number of Cells = ' num2str(length(cells2useAbs))])
    zpeaks2 = zpeaks2(cells2use);
    bistratType = bistratType(cells2use);

    %---------- load z-distributions ----------%
    zProfile = zeros(length(cells2use), 120);
    cellNames = cell(length(cells2use),1);
    n = 0;
    for i = 1:length(cells2useAbs)
        % find cell
        ret = retina{1,2}(cells2useAbs(i));
        if ret == 1
            ret = 'L';
        else
            ret = 'R';
        end
        num = ['000',int2str(cells2useAbs(i))]; num = num(end-3:end);
        exp = ['0000',int2str(experiment{1,2}(cells2useAbs(i)))]; exp = exp(end-4:end);
        sc = ['0',int2str(scan{1,2}(cells2useAbs(i)))]; sc = sc(end-1:end);
        ce = ['0',int2str(cellNumber{1,2}(cells2useAbs(i)))]; ce = ce(end-1:end);
        cellName = [num,'_',exp,'_',int2str(mouse{1,2}(cells2useAbs(i))),ret,'_C',sc,'_',ce];
        cellNames{i} = cellName;

        clear arborDensity
        eval(['clear ' arbourName])
        load(fullfile(dname ,[cellName,'_zDist']), arbourName);
        
            disp(int2str(cells2useAbs(i)))
            
        eval(['arborDensity = ' arbourName ';']);
        
        n = n + 1;
        if i>1 && length(arborDensity(:))~=size(aD,2)
            aD(n,:)=zeros(1,size(aD,2));
        else
            aD(n,:) = arborDensity(:)*(sqrt(areas(cells2use(i)))/umPerPix)/norm(arborDensity(:));
        end
        zProfile(n,:) = sum(sum(arborDensity,2),3)/sum(arborDensity(:));
    end

    %---------- Clean Arbor Densities based on zProfiles ----------%
    OzProfile = zProfile./repmat(max(zProfile')',1,size(zProfile,2));
    X = -30:.5:29.5; X = X./12;
    for i = 1:size(zProfile,1)
         %----- Smooth -----%           
         zProfile(i,:) = movmean(OzProfile(i,:),5);                  
         
         %----- Threshold -----%
%          T1 = .5; %.3679;
         if bistratType(i) == 1
            T2 = zProfile(i,find(X > zpeaks2(i),1,'first')) * 0.5;
            if T2 < T1; TT = T2; else; TT = T1; end                            
            idx = find(zProfile(i,:) < TT); 
         else
            idx = find(zProfile(i,:) < T1);
         end                  
         zProfile(i,idx) = 0;
         dp = reshape(aD(i,:),120,xSize,xSize);         
         dp(idx,:,:) = 0;
                  
         %----- DownSample in Z -----%
         if downZ == 1              
             downsamplerate = 120/nz2;
             dp2 = zeros(nz2, xSize, xSize);
             for j = 1:xSize
                 for k = 1:xSize
                    dp2(:, j, k) = resample(dp(:,j,k), nz2 ,120);
                 end
             end
             dp = dp2;
             clear dp2
         else
             nz2 = 120;
         end         
         npoints = nz2 * xSize * xSize;
                    
             
         %----- Invert Bistratified -----%
         if i == 1
             aD2 = zeros(length(cells2useAbs), npoints);
         end
         if invertBistrat == 1 && bistratType(i) == 1
             aD2(i,:) = -reshape(dp,1,npoints); 
         else
            aD2(i,:) = reshape(dp,1,npoints);    
         end
         
         
         
         %----- Visualize ------%
%          if bistratType(i) == 1
%             plot(X,mean(mean(dp,3),2)./max(dp(:)),'color', [0 .447 .741],'linewidth',2); hold on
%             plot(X,OzProfile(i,:),'color', [.85 .325 .098]);
%             line(zeros(101,1),0:.01:1); line(ones(101,1),0:.01:1); 
%             title(['Cell: ' num2str(i) '.   TwoPeaks: ' num2str((zpeaks2(i)))]); hold off
%             reply = input('Do you want more? Y/N [Y]:','s');
%          end
%             
    end
    cData = aD2;
    clear  aD2

%% Dimenstionality Reduction 
% PCA and tSNE see ok. sPCA and NNMF seem nonsense. 5 - 25 dimensions.
if doDimReduction == 1
    
    switch dimreductionMeth
        %----- Do PCA -----% 
        case 'PCA'
            [dData, ~] = pca(cData', 'Algorithm', 'svd' ,'NumComponents', NumberDimensions);
            
        %----- Do SVD -----% 
        case 'SVD'
            
        %----- Do Sparse PCA -----%    
        case 'sPCA' 
        [dData ~] = spca(cData', [], NumberDimensions, inf, 1,  1000, 1e-6, true);
    
        %----- Do NNMF -----% 
        case 'NNMF'
            [dData, ~] = nnmf(cData, NumberDimensions);    
            
        %----- Do tSNE -----%
        case 'tSNE'
            dData = tsne(cData,'Algorithm','exact','Distance','euclidean', 'NumDimensions', NumberDimensions);
            
        otherwise
            
            
    end    
    Data2Cluster = dData;
else
    Data2Cluster = cData;
    
end
% size(Data2Cluster)
% plotmatrix(Data2Cluster)

    
%% Cluster

disp('     E - Linkage ...'); 
%K = 15;
alpha = 1;
%SimIndex = zeros(30,size(aD,1));
%clusterLUT = zeros(30,size(aD,1));
Dmat = pdist(aD,'euclidean');
Link = elinkage(Dmat, alpha);

d = size(aD,1);
clusterLUT = zeros(30, d);
SimIndex =  zeros(30, d); %put in decision matrix
RanIndex = zeros(30,d);

for K = 1:30
%K = 3;
originalCluster = cluster(Link, 'maxclust',K);
%BB=[94:99 142 153 170 186 201 209 240 246 259 264 285 290 316 318 321 328 354];
tic


parfor i = 1:size(aD,1)
temp = aD;
knownCell = aD(i,:);    
temp(i,:) = []; %%take out the cell from the original dataset to recluster
knownCellId = originalCluster(i);
oldCluster = originalCluster;
oldCluster(i) = [];  %take out the test cell 
%oldCluster = [oldCluster; knownCellId]; %%move the test cell to bottom
knownCellPos = i;
%%%recluster the remaining 362 cells
tempDmat = pdist(temp, 'euclidean');
myLinkage = elinkage(tempDmat, alpha);
%%%%%%%%%%%%%Cutoff, check for similarity and RI%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxCutLevel = 1;
% normalize linkage values
myLinkage(:,3) = myLinkage(:,3)/myLinkage(end,3);
% cutting levels for the linkage
cutStep = 0.001; threshold = cutStep:cutStep:maxCutLevel;
% results of all cuts
randIndices = zeros(size(threshold)); adjustedRandIndices = zeros(size(threshold)); typeConfs = zeros(size(threshold)); widths = zeros(size(threshold));
rowSplits = zeros(size(threshold)); columnSplits = zeros(size(threshold));
% initialize hard-clustering variable for all cutting levels (allT)
allT=ones(size(myLinkage,1)+1,numel(threshold));

for kk = 1:numel(threshold)
    % cut the dendrogram at threshold(kk)
    T = cluster(myLinkage,'cutoff',threshold(kk),'criterion','distance');
%     distance = zeros(max(T),1);
%     %[S,I] = sort(T,'ascend');
%     for cnum = 1:max(T)
%         ClusterIndex = find(T==cnum);
%         Y = size(ClusterIndex,1);
%         TempCluster = zeros(Y, size(temp,2));
%         
%         for y = 1:Y
%             TempCluster(y,:) = temp(ClusterIndex(y),:); %
%         end
%         if Y == 1
%             ADMean = TempCluster;
%         else
%             ADMean = mean(TempCluster);
%         end
%         distance(cnum) = pdist2(knownCell, ADMean, 'euclidean'); 
%     end
%     [mindis, groupCellBelongTo] = min(distance);
%     Tnew = [T; groupCellBelongTo];
    % calculate clustering accuracy metrics
   allT(:,kk) = T; 
   [AR,RI,MI,HI,rS,cS,typeConfusions]=reportConfusionsAndRI(T,oldCluster);
   randIndices(kk) = RI; adjustedRandIndices(kk) = AR; rowSplits(kk) = rS; columnSplits(kk) = cS; typeConfs(kk) = typeConfusions;
end
% calculate the width of the range of cutting levels yielding the exact same clustering, for each cutting level
for kk2 = 1:numel(threshold)
    flag=true;width=1;while flag&&(kk2-width>0);[AR,RI1,MI,HI,~,~,~]=reportConfusionsAndRI(allT(:,kk2-width),allT(:,kk2));
        if RI1==1;width=width+1;else;flag=false;end;end;
    flag=true;tmpwidth=0;while flag&&(kk2+tmpwidth+1<=numel(threshold));[AR,RI2,MI,HI,~,~,~]=reportConfusionsAndRI(allT(:,kk2+tmpwidth+1),allT(:,kk2));
        if RI2==1;tmpwidth=tmpwidth+1;else;flag=false;end;end;
    widths(kk2) = width+tmpwidth;
end

[mini0, pos0] = min(typeConfs); allMin = find(typeConfs==mini0); % among minimum type-confusions ...
[~, postemp] = max(randIndices(allMin)); pos0=allMin(postemp); th0=threshold(pos0); width=widths(pos0);
myRowSplits = rowSplits(pos0); myColumnSplits = columnSplits(pos0);
T = cluster(myLinkage,'cutoff',th0,'criterion','distance');
distance = zeros(max(T),1);

for cnum = 1:max(T)
    ClusterIndex = find(T==cnum);
    Y = size(ClusterIndex,1);
    TempCluster = zeros(Y, size(temp,2));
    %
    for y = 1:Y
        TempCluster(y,:) = temp(ClusterIndex(y),:); %
    end
    if Y == 1
        ADMean = TempCluster;
    else
        ADMean = mean(TempCluster);
    end
    distance(cnum) = pdist2(knownCell, ADMean, 'euclidean');
end
[mindis, groupCellBelongTo] = min(distance);

NewClusterGroupIndex = find(T==groupCellBelongTo);
OldClusterGroupIndex = find(oldCluster==knownCellId); 
%     Tnew = [T; groupCellBelongTo];

JaccardIndex = jaccard(NewClusterGroupIndex, OldClusterGroupIndex);
RanIndex(K,i) = randIndices(pos0);
SimIndex(K,i) =  JaccardIndex; %put in decision matrix
clusterLUT(K,i) = numel(unique(T));
end
toc
end


%%%post process
% option = [];
% count = 1;
% for i = 1:30
%     groupInTest = clusterLUT(i,1);
%     if size(find(clusterLUT(i,:) ~= groupInTest),2) == 0
%         option(count) = groupInTest;
%         count = count + 1;
%     end
% end
% %%%%%%%choose cluster
% display(option)

Score = [];
for i = 1:30
    res = zeros(30,size(clusterLUT,2));
    [r,c]= find(clusterLUT==i);
    res(r,c) = 1;
    resFin = res.* RanIndex;
    Score(i) = sum(resFin(:));
end

Matrix = [];
Score = round(Score);
for i = 1:30
    for j = 1:Score(i)
        Matrix = [Matrix;i];
    end
end


save("OurDataClusterQDScript_02.mat", "clusterLUT", "Matrix", "SimIndex", "RanIndex");


histogram(Matrix)
yt = get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',round(100*yt/size(Matrix,1)))
xlabel("Number of Clusters")
ylabel("Percentage")

%% Validate Clustering: Pick Number of Clusters

%---------- Validation Indexes ----------%

% save(fullfile(path2save,[name2save,'_n',int2str(normalize) ]),'C','CMeths','K','OptimalK','aD','cells2useAbs','cells2use','params','arbourName','normalize','VI','AR','RI','MI','HI')
end