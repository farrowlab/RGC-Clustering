function     CheckClassificationRGCarborDensitiesQD(xSize,invertBistrat, includeMaybes,types2analyse,dname, dbase,path2save,downZ,nz2,doDimReduction,dimreductionMeth,NumberDimensions,CMeths,dist,arbourName,name2save,zProfileThreshold)
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

    
%% Check Classify
%labels = load('/home/quan/Desktop/ClassificationScripts/classesLABEL6.txt');
labels = load('/home/quan/Desktop/ClassificationScripts/NewClassificationParams/classesLABEL.txt')
count = 0;


zProfile = zProfile./repmat(max(zProfile')',1,size(zProfile,2));
sProfile = zProfile';

% % %%%%Remove the one that has a prob lower than 30%, the noisy one
 % indexToDelete = [];
 % for i = 1:size(labelsProb,1)
 %     if labels(i) ~= labels2(i)
% %  %if max(labelsProb(i,:)) < 0.2 
%          indexToDelete = [indexToDelete, i];
% %   %      %sProfile(:,i) = []; 
% %         %labels(i) = [];
% %  %   elseif (max(labelsProb2(i,:)) + max(labelsProb(i,:)))/2 < 0.2
% %  %       indexToDelete = [indexToDelete, i];    end
%      end
%  end
%  sProfile(:, indexToDelete) = [];
%labels(indexToDelete) = 0;
%fileName(indexToDelete);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avProfile = zeros(120,15);
for i = 1:15
    tmp = find(labels==i);
    data = sProfile(:,tmp);
    data = data./repmat(max(data),120,1);
    av = mean(data');
    avProfile(:,i)=av;
end
avProfile=avProfile./repmat(max(avProfile),120,1);

figure
X = -30:.5:29.5; X = X./12;
cmap = colormap(lines(15));
step = 1.1 * max(avProfile(:));
hold on
for k = 1:15
    plot(X, step * (k - 1) + ((avProfile(:,k)')), 'k'); hold on
    
    text(-2.7,step * (k - 1) ,int2str(k))
    stepmax = step + step
end
xlim([-2 2.4])
set(gca,'Box','off', 'YColor','w')
ylim([0 15 * step])

end