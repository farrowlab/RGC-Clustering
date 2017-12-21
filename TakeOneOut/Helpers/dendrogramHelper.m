%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function labelDendrogram(thisLinkage,labels,colorLabels,saveName)
figure;h = dendrogram(thisLinkage,0);
% cell IDs in the dendrogram leaves
curLabels=get(gca,'XTickLabel'); for kk = 1:size(curLabels,1); numCurLabels(kk) = str2num(curLabels(kk,:)); end;
% label leaves
for kk = 1:size(curLabels,1);
  x = kk; y = -0.01; t=text(x,y,labels{numCurLabels(kk)},'Color',colorLabels{numCurLabels(kk)});
  set(t,'HorizontalAlignment','right','VerticalAlignment','top'); %,'Rotation',90);
end
%Prettier
set(gca,'XLim',[0,size(thisLinkage,1)+2],'YLim',[0,1.05], 'XTickLabel', {}, 'XTick',[], 'YTickLabel', {'0','0.2','0.4','0.6','0.8','1'}, 'YTick',[0:0.2:1]);
for kk = 1:numel(h); set(h(kk),'Color','black','LineWidth',2); end;
pixelsPerInch = 70; set(gcf,'Color','w','PaperUnits', 'inches', 'PaperPosition',[100 100 4000 800]/pixelsPerInch); set(gca,'box','off');
print(gcf,'-dpng',sprintf('-r%d',pixelsPerInch), strcat(saveName,'_dendrogram.png')); close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function labels=readLabels
% cell types of the 363 cells of the dataset as presented in the paper
labels=cell(364,1); labels{364}='N';
colorLabels=cell(364,1); colorLabels{364}='r';
AA=[84:93 127 304 355]; BB=[94:99 142 153 170 186 201 209 240 246 259 264 285 290 316 318 321 328 354];
CC=[100:105 125 126 130 136 137 149 150 151 172 173 195 205 206 214 215 228 235 237 238 293 299 305 315 319 342 344 345 346 347 348 356];
DD=[65:83 112 119 131 133 145 154 156 157 162 163 164 168 169 171 177 178 181 185 190 191 194 197 198 202 210 218];
DD=[DD 220 221 239 241 245 247 250 251 253 255 268 272 274 286 287 296 306 307 309 310 327 329 343 350 352];
EE=[42:57 108:111 216]; FF=[33:41 114 115 159 160 176 203 225 242 265 276 300 311 330];
GG=[58:64 263 337 363];
HH=[1:32 200 323];
II=[106:107 117 121 122 123 128 141 158 166 179 219 222 226 231 243 244 266 269 288 292 298];
VV=[193 224 283 336 358 359 362];
UU=[129 146 161 167 175 187 189 192 196 204 212 213 229 232 254 281 282 308 312 314 320 322 331 334 349 360];
YY=[139 140 147 148 174 217 227 230 234 252 267 273 277 289 294 295 297 340];
WW=[258 270 280 303 325 341 353 357];
XX=[113 116 132 134 138 152 155 165 180 182 183 184 188 207 233 236 248 256 257 271 275 278 284 302 313 324 326 333 338 339];
ZZ=[118 120 124 135 143 144 199 208 211 223 249 260 261 262 279 291 301 317 332 335 351 361];
labels(AA)={'a'}; labels(BB)={'b'}; labels(CC)={'c'}; labels(DD)={'d'}; labels(EE)={'e'}; labels(FF)={'f'}; labels(GG)={'g'}; labels(HH)={'h'};
labels(II)={'i'}; labels(UU)={'u'}; labels(VV)={'v'}; labels(WW)={'w'}; labels(XX)={'x'}; labels(YY)={'y'}; labels(ZZ)={'z'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hardClusterStats = cutDendrogram(myLinkage,clusterIDsForKnownCells,knownCellPositions,maxCutLevel)
% clusterIDsForKnownCells: array of integers of length numel(knowmCellPositions), where cells with same IDs belong to the same cluster
if nargin < 4; maxCutLevel = 0.1; end;
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
 allT(:,kk) = T; [AR,RI,MI,HI,rS,cS,typeConfusions]=reportConfusionsAndRI(T(knownCellPositions),clusterIDsForKnownCells);
 randIndices(kk) = RI; adjustedRandIndices(kk) = AR; rowSplits(kk) = rS; columnSplits(kk) = cS; typeConfs(kk) = typeConfusions;
end
% calculate the width of the range of cutting levels yielding the exact same clustering, for each cutting level
for kk = 1:numel(threshold)
 flag=true;width=1;while flag&&(kk-width>0);[AR,RI1,MI,HI,~,~,~]=reportConfusionsAndRI(allT(:,kk-width),allT(:,kk));
 if RI1==1;width=width+1;else;flag=false;end;end;
 flag=true;tmpwidth=0;while flag&&(kk+tmpwidth+1<=numel(threshold));[AR,RI2,MI,HI,~,~,~]=reportConfusionsAndRI(allT(:,kk+tmpwidth+1),allT(:,kk));
 if RI2==1;tmpwidth=tmpwidth+1;else;flag=false;end;end;
 widths(kk) = width+tmpwidth;
end

[mini0, pos0] = min(typeConfs); allMin = find(typeConfs==mini0); % among minimum type-confusions ...
[~, postemp] = max(randIndices(allMin)); pos0=allMin(postemp); th0=threshold(pos0); width=widths(pos0);
myRowSplits = rowSplits(pos0); myColumnSplits = columnSplits(pos0);
T = cluster(myLinkage,'cutoff',th0,'criterion','distance');

hardClusterStats.minTypeConfs = mini0;
hardClusterStats.structuralSplits = myRowSplits;
hardClusterStats.geneticSplits = myColumnSplits;
hardClusterStats.minCutForMinTypeConfs = th0;
hardClusterStats.cutWidthForMinTypeConfs = width*cutStep;
hardClusterStats.ClusterCountAtMinCut = numel(unique(T));
hardClusterStats.hardClusters = T;
hardClusterStats.riAtMinCut = randIndices(pos0);
hardClusterStats.ariAtMinCut = adjustedRandIndices(pos0);
hardClusterStats.ri = randIndices;
hardClusterStats.ari = adjustedRandIndices;
hardClusterStats.allTypeConfs = typeConfs;
hardClusterStats.allCutWidths = widths*cutStep;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AR,RI,MI,HI,rowSplits,columnSplits,typeConfusions] = reportConfusionsAndRI(c1,c2)

% This code is slightly modified from the original by Uygar Sümbül.
% The copyright information for the original file is generated below.
% The original file can be downloaded from the following website as of Dec. 23, 2013:
% http://www.mathworks.com/matlabcentral/fileexchange/13916-simple-tool-for-estimating-the-number-of-clusters/content/valid_RandIndex.m

%RANDINDEX - calculates Rand Indices to compare two partitions
% ARI=RANDINDEX(c1,c2), where c1,c2 are vectors listing the
% class membership, returns the "Hubert & Arabie adjusted Rand index".
% [AR,RI,MI,HI]=RANDINDEX(c1,c2) returns the adjusted Rand index,
% the unadjusted Rand index, "Mirkin's" index and "Hubert's" index.
%
% See L. Hubert and P. Arabie (1985) "Comparing Partitions" Journal of
% Classification 2:193-218

%(C) David Corney (2000) D.Corney@cs.ucl.ac.uk

if nargin < 2 | min(size(c1)) > 1 | min(size(c2)) > 1
   error('RandIndex: Requires two vector arguments')
   return
end

C=Contingency(c1,c2); %form contingency matrix

n=sum(sum(C));
nis=sum(sum(C,2).^2); %sum of squares of sums of rows
njs=sum(sum(C,1).^2); %sum of squares of sums of columns

t1=nchoosek(n,2); %total number of pairs of entities
t2=sum(sum(C.^2)); %sum over rows & columnns of nij^2
t3=.5*(nis+njs);

%Expected index (for adjustment)
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3; %no. agreements
D= -t2+t3; %no. disagreements

if t1==nc
   AR=0; %avoid division by zero; if k=1, define Rand = 0
else
   AR=(A-nc)/(t1-nc); %adjusted Rand - Hubert & Arabie 1985
end

RI=A/t1; %Rand 1971 %Probability of agreement
MI=D/t1; %Mirkin 1970 %p(disagreement)
HI=(A-D)/t1; %Hubert 1977 %p(agree)-p(disagree)

tmp1=sum(C>0,1); tmp2=sum(C>0,2);
rowSplits = sum(tmp1)-nnz(tmp1); columnSplits = sum(tmp2)-nnz(tmp2);
typeConfusions = rowSplits+columnSplits; % number of splits in both directions!
%for kk=1:size(C,2)
% [maxi,pos]=max(C(:,kk)); if pos>kk; tmp = C(pos,:); C(pos,:) = C(kk,:); C(kk,:) = tmp; end;
%end
%tmp = diag(diag(C)); tmp = [tmp; zeros(size(C,1)-size(tmp,1),size(tmp,2))]; tmp = [tmp zeros(size(tmp,1), size(C,2)-size(tmp,2))]; tmp = C-tmp;
%typeConfusions = nnz(tmp);

function Cont=Contingency(Mem1,Mem2)

if nargin < 2 | min(size(Mem1)) > 1 | min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end
