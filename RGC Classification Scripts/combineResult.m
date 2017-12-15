function combineResult
clear all
labels = load('/media/areca_raid/Classification Scripts/Labels/classesLABEL.txt');
labelsProb = load('/media/areca_raid/Classification Scripts/Labels/classesPROB.txt');

labels2 = load('/media/areca_raid/Classification Scripts/Labels/classesLABEL2.txt');
labelsProb2 = load('/media/areca_raid/Classification Scripts/Labels/classesPROB2.txt');

%labels3 = load('/home/quan/Desktop/ClassificationScripts/classesLABEL3.txt');
%labelsProb3 = load('/home/quan/Desktop/ClassificationScripts/classesPROB3.txt');

%labels4 = load('/home/quan/Desktop/ClassificationScripts/RF,cglptuwx6classesLABEL.txt');
%labelsProb4 = load('/home/quan/Desktop/ClassificationScripts/RFclassesPROB.txt');

files = dir('/media/areca_raid/LabPapers/SCRouter/Data/dendriticArbors/*_zDist.mat');
cd('/media/areca_raid/LabPapers/SCRouter/Data/dendriticArbors/');

n = 1;
zProfile = [];

%fileName = cell(1,325);
%%%%%%%%%%resize and create zProfile%%%%%%%%%%%%%%%
for file = files'
   % if count == 325
    %    break
   % else
   %     count = count+1;
   % end
    %resultname = [resultname;file.name];
    a = load(file.name,'arborDensity2');
    a = a.arborDensity2;
    %fileName{count}= file.name;
    a = permute(a, [3,2,1]);
    %  [x, y, z] = size(a);
    %  X = [X;x];
    %  Y = [Y;y];
    %  Z = [Z;z];
    ny=29;nx=29;nz=120; %% desired output dimensions
    [y x z]=...
        ndgrid(linspace(1,size(a,1),ny),...
        linspace(1,size(a,2),nx),...
        linspace(1,size(a,3),nz));
    aOut=interp3(double(a),x,y,z);
    
    a = permute(aOut, [3,2,1]);
    
    zProfile(:,n) = sum(sum(a,2),3)/sum(a(:));
    n = n+1;

end
zProfile = zProfile./repmat(max(zProfile')',1,size(zProfile,2));
sProfile = zProfile;

% % %%%%Remove the one that has a prob lower than 30%, the noisy one
  indexToDelete = [];
  for i = 1:size(labelsProb,1)
      if labels(i) ~= labels2(i)
% %  %if max(labelsProb(i,:)) < 0.2 
          indexToDelete = [indexToDelete, i];
% %   %      %sProfile(:,i) = []; 
% %         %labels(i) = [];
% %  %   elseif (max(labelsProb2(i,:)) + max(labelsProb(i,:)))/2 < 0.2
% %  %       indexToDelete = [indexToDelete, i];    end
      end
  end
%  sProfile(:, indexToDelete) = [];
labels(indexToDelete) = 0;
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
