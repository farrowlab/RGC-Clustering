function createTestDataset
files = dir('/media/areca_raid/LabPapers/SCRouter/Data/dendriticArbors/*_zDist.mat');
cd('/media/areca_raid/LabPapers/SCRouter/Data/dendriticArbors/');
ADtestresult = [];
resultname = [];
n = 1;

for file = files'
    %resultname = [resultname;file.name];
    a = load(file.name,'arborDensity2');
    a = a.arborDensity2;
    [r,c,h] = size(a);
    
    a = permute(a, [3,2,1]);
    ny=29;nx=29;nz=120; %% desired output dimensions
    [y x z]=...
        ndgrid(linspace(1,size(a,1),ny),...
        linspace(1,size(a,2),nx),...
        linspace(1,size(a,3),nz));
    aOut=interp3(double(a),x,y,z);
    
    
    a = permute(aOut, [3,2,1]);
    a = mat2gray(a);
    
    zProfile(:,n) = sum(sum(a,2),3)/sum(a(:));
    n = n+1;
    
    aOut = a(:)';
    area = [area; trapz(aOut)];
    %ADtestresult = [ADtestresult; aOut];
end
%[testdData,~] = pca(ADtestresult', 'Algorithm','svd','NumComponents', 1000);
zProfile = zProfile./repmat(max(zProfile')',1,size(zProfile,2));

for i = 1:size(zProfile,2)
    zMean_test(i,:) = mean(zProfile(:,i));
    zVariance_test(i,:) = var(zProfile(:,i));
    [peak, loc, w, p] = findpeaks(zProfile(:,i));
    [mpeak, midx] = max(peak);
    mwidth = w(midx);
    zWidth_test1(i,:) = mwidth;
    zLocation_test1(i,:) = loc(midx);
    
    peak(midx) = []; %% delete that max
    w(midx) = [];
    loc(midx) = [];
    
    %%now second max
    [mpeak2, midx2] = max(peak); %%get second max peak
    mwidth2 = w(midx2);
    zWidth_test2(i,:) = mwidth2;
    zLocation_test2(i,:) = loc(midx2);
    
    
    peak(midx2) = [];
    w(midx2) = [];
    loc(midx2) = [];
    
    %third max
    [mpeak3, midx3] = max(peak); %%get second max peak
    mwidth3 = w(midx3);
    zWidth_test3(i,:) = mwidth3;
    zLocation_test3(i,:) = loc(midx3);
    
    peak(midx3) = [];
    w(midx3) = [];
    loc(midx3) = [];
    
    %third max
    [mpeak4, midx4] = max(peak); %%get second max peak
    mwidth4 = w(midx4);
    zWidth_test4(i,:) = mwidth4;
    zLocation_test4(i,:) = loc(midx4);
    
    peak(midx4) = [];
    w(midx4) = [];
    loc(midx4) = [];
    
    %third max
    [mpeak5, midx5] = max(peak); %%get second max peak
    mwidth5 = w(midx5);
    zWidth_test5(i,:) = mwidth5;
    zLocation_test5(i,:) = loc(midx5);
end
testdata = [zLocation_test1, zWidth_test1, zLocation_test2, zWidth_test2,zLocation_test3, zWidth_test3, zLocation_test4, zWidth_test4,zLocation_test5, zWidth_test5, zMean_test, zVariance_test, area];
cd('/media/areca_raid/Classification Scripts/Data');
dlmwrite('testdata.txt',testdata);
end
