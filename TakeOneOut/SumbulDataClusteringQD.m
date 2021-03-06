clear
load('/media/areca_raid/LabPapers/SCRouter/rgcdata/Sumbul/arborDensities.mat')

aD = arborDensities';


disp('     E - Linkage ...'); 
%K = 15;
alpha = 1;

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

    % calculate clustering accuracy metrics
   allT(:,kk) = T; 
   [AR,RI,MI,HI,rS,cS,typeConfusions]=reportConfusionsAndRI(T,oldCluster);
   randIndices(kk) = RI; adjustedRandIndices(kk) = AR; rowSplits(kk) = rS; columnSplits(kk) = cS; typeConfs(kk) = typeConfusions;
end
% calculate the width of the range of cutting levels yielding the exact same clustering,for each cutting level
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



% %%%post process
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


save("SumbulDataClusterQD.mat", "clusterLUT", "Matrix", "SimIndex", "RanIndex");


histogram(Matrix)
yt = get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',round(100*yt/size(Matrix,1)))
xlabel("Number of Clusters")
ylabel("Percentage")






