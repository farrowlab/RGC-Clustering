function     classifyRGCarborDensitiesQD(xSize,invertBistrat, includeMaybes,types2analyse,dname, dbase,path2save,downZ,nz2,doDimReduction,dimreductionMeth,NumberDimensions,CMeths,dist,arbourName,name2save,zProfileThreshold)
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

    
%% Classify

disp('     E - Linkage ...'); 
%K = 15;
alpha = 1;
%SimIndex = zeros(30,size(aD,1));
%clusterLUT = zeros(30,size(aD,1));
[a, b] = size(aD);
ADresult = [];
zProfile2 = [];
n=1;
area2 = [];
for i = 1:a
    temp = reshape(aD(i,:),[120, 28 ,28]);
    
    temp = permute(temp, [3,2,1]);
    
    ny=29;nx=29;nz=120; %% desired output dimensions
    [y x z]=...
        ndgrid(linspace(1,size(temp,1),ny),...
        linspace(1,size(temp,2),nx),...
        linspace(1,size(temp,3),nz));
    aOut=interp3(double(temp),x,y,z);
    temp = permute(aOut, [3,2,1]);
    temp = mat2gray(temp);
    zProfile2(:,n) = sum(sum(temp,2),3)/sum(temp(:));
    n = n+1;
   aOut = temp(:)';
   area2 = [area2; trapz(aOut)];
  % ADresult = [ADresult; aOut];
end
zProfile2 = zProfile2./repmat(max(zProfile2')',1,size(zProfile2,2));
%[dData,~] = pca(ADresult', 'Algorithm','svd','NumComponents', 25);

for i = 1:size(zProfile2,2)
     zMean_train(i,:) = mean(zProfile2(:,i)); %mean
     zVariance_train(i,:) = var(zProfile2(:,i)); %variance
     [peak, loc, w, p] = findpeaks(zProfile2(:,i)); %peak location
     
     %% first max 
     [mpeak, midx] = max(peak);
     mwidth = w(midx);
     zWidth_train1(i,:) = mwidth;
     zLocation_train1(i,:) = loc(midx);
     
     peak(midx) = []; %% delete that max
     w(midx) = [];
     loc(midx) = [];
     
     %% second max 
     [mpeak2, midx2] = max(peak); %%get second max peak
     mwidth2 = w(midx2);
     zWidth_train2(i,:) = mwidth2;
     zLocation_train2(i,:) = loc(midx2);
     
     peak(midx2) = [];
     w(midx2) = [];
     loc(midx2) = [];
     
     %%third max
     [mpeak3, midx3] = max(peak);
     mwidth3 = w(midx3);
     zWidth_train3(i,:) = mwidth3;
     zLocation_train3(i,:) = loc(midx3);
     
     peak(midx3) = [];
     w(midx3) = [];
     loc(midx3) = [];
     
     %%fourth max
     [mpeak4, midx4] = max(peak);
     mwidth4 = w(midx4);
     zWidth_train4(i,:) = mwidth4;
     zLocation_train4(i,:) = loc(midx4);
     
     peak(midx4) = [];
     w(midx4) = [];
     loc(midx4) = [];
     
     %%fifth max
     [mpeak5, midx5] = max(peak);
     mwidth5 = w(midx5);
     zWidth_train5(i,:) = mwidth5;
     zLocation_train5(i,:) = loc(midx5);
end

path1='/media/areca_raid/LabPapers/SCRouter/Katja';

cd(path1)


cd('/media/areca_raid/RGC-Clustering/RGC Classification Scripts/');
traindata = [zLocation_train1, zWidth_train1, zLocation_train2, zWidth_train2, zLocation_train3, zWidth_train3, zLocation_train4, zWidth_train4,zLocation_train5, zWidth_train5,zMean_train, zVariance_train,area2];
dlmwrite('/media/areca_raid/RGC-Clustering/RGC Classification Scripts/Results/OurDataTest.txt',traindata);   


% option = [];
% count = 1;
% for i = 1:30
%     groupInTest = clusterLUT(i,1);
%     if size(find(clusterLUT(i,:) ~= groupInTest),2) == 0
%         option(count) = groupInTest;
%         count = count + 1;
%     end
% end

%% Validate Clustering: Pick Number of Clusters

%---------- Validation Indexes ----------%

% save(fullfile(path2save,[name2save,'_n',int2str(normalize) ]),'C','CMeths','K','OptimalK','aD','cells2useAbs','cells2use','params','arbourName','normalize','VI','AR','RI','MI','HI')
end
