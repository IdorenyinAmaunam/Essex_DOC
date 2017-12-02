function [r p] = findCorrWithFMAdiff(Map,FMAdiff, SubID, Side, IndBCI, IndSham)

SMap = nan(size(Map{1},2),size(Map{1},1),0);
for sub=1:length(SubID)

    % Mirror as needed
    if(strcmp(Side{sub},'R'))
        doflip = 1;
    else
        doflip = 0;
    end
    
    Map{sub} = flipSide16Matrix(Map{sub}',doflip);
    
    % Find average-based metrics
    AvgBestERD = min(Map{sub}(:));
    AvgBestERS = max(Map{sub}(:));
    AvgBestAbs = max(abs(Map{sub}(:)));
    Sorted = sort(Map{sub}(:),'ascend');
    SortedAbs = sort(abs(Map{sub}(:)),'descend');
    AvgBest5ERD = mean(Sorted(1:5));
    AvgBest5ERS = mean(Sorted(end-4:end));
    AvgBest5Abs = mean(SortedAbs(1:5));
    AvgContraMu = matmean(Map{sub}([2 3 7 8 12 13],[3:6]));
    AvgContraBeta = matmean(Map{sub}([2 3 7 8 12 13],[8:14]));
    AvgIpsiMu = matmean(Map{sub}([5 6 10 11 15 16],[3:6]));
    AvgIpsiBeta = matmean(Map{sub}([5 6 10 11 15 16],[8:14]));
    AvgMedMu = matmean(Map{sub}([4 9 14],[3:6]));
    AvgMedBeta = matmean(Map{sub}([4 9 14],[8:14]));    
       
    Metrics(sub,:) = [AvgBestERD AvgBestERS AvgBestAbs ...
        AvgBest5ERD AvgBest5ERS AvgBest5Abs ...
        AvgContraMu AvgContraBeta ...
        AvgIpsiMu AvgIpsiBeta ...
        AvgMedMu AvgMedBeta];    
    
    SMap = cat(3,SMap,Map{sub});
end


%% Nat Comm revision reply 1.13
% Colors for BCI and Sham
Color = {[255 40 40]/255,[102 184 184]/255};
Metric = abs(Metrics(:,8));

close all;
figure(103);
bar(1,mean(Metric(IndBCI)),'FaceColor',Color{1});hold on;
bar(2,mean(Metric(IndSham)),'FaceColor',Color{2});
errorbar(1, mean(Metric(IndBCI)), std2(Metric(IndBCI)),'.k','LineWidth',3);
errorbar(2, mean(Metric(IndSham)), std2(Metric(IndSham)),'.k','LineWidth',3);
[h p] = ttest2(Metric(IndBCI),Metric(IndSham));
axis([0 3 -0.05 0.8]);
h = sigstar({[1,2]},[p]);
set(h,'LineWidth',3);
hold off;
legend({'BCI','Sham'},'FontSize',20);
set(gca,'XTick',[1:2]);
set(gca,'XTickLabel',{'BCI' , 'Sham'},'FontSize',20);
ylabel('Fugl-Meyer Assessment Score (FMA) POST-PRE','FontSize',20);
xlabel('|ERSP_{Hit} - ERSP_{Miss}|, contralesional \beta band [18,24] Hz)','FontSize',20);
set(gca,'FontSize',20);
set(gca,'LineWidth',3);



figure(104);
h = plot(Metric,FMAdiff,'.','Color','k','MarkerSize',1);
hold on;
h1 = plot(Metric(IndBCI),FMAdiff(IndBCI)','.','Color',Color{1},'MarkerSize',30);
h2 = plot(Metric(IndSham),FMAdiff(IndSham)','.','Color',Color{2},'MarkerSize',30);
hls = lsline;
set(hls,'LineWidth',3);
hold off;
axis([-0.05 1.05 -3 20]);
legend([h1 h2],{'BCI','Sham'},'Location','NorthWest','FontSize',20);
ylabel('Fugl-Meyer Assessment Score (FMA) POST-PRE','FontSize',20);
xlabel('|ERSP_{Hit} - ERSP_{Miss}|, contralesional \beta band [18,24] Hz)','FontSize',20);
set(gca,'FontSize',20);
[r p] = corr(Metric,FMAdiff')



% r = [];
% p = [];
% for m=1:size(Metrics,2)
%     [r(m) p(m)] = corr(Metrics(:,m),FMAdiff');
% end
% 
% rr=[];
% pp=[];
% for i=1:size(SMap,1)
%     for j=1:size(SMap,2)
%         [rr(i,j) pp(i,j)] = corr(squeeze(SMap(i,j,:)),FMAdiff');        
%     end
% end
% 
% rrr=rr;
% rrr(pp>=0.05)=0;
% figure(1);
% imagesc(rrr);
% 
% i=7; % 8
% j=8; % 12
% figure(2);
% %plot(squeeze(SMap(i,j,:)),FMAdiff','*k',squeeze(SMap(i,j,IndBCI)),FMAdiff(IndBCI)','*b',squeeze(SMap(i,j,IndSham)),FMAdiff(IndSham)','*r');
% plot(abs(squeeze(SMap(i,j,:))),FMAdiff','*k',abs(squeeze(SMap(i,j,IndBCI))),FMAdiff(IndBCI)','*b',abs(squeeze(SMap(i,j,IndSham))),FMAdiff(IndSham)','*r');
% lsline;
% 
% k=8;
% figure(3);
% plot(Metrics(:,k),FMAdiff','*k',Metrics(IndBCI,k),FMAdiff(IndBCI)','*b',Metrics(IndSham,k),FMAdiff(IndSham)','*r');
% lsline;
% figure(4);
% plot(abs(Metrics(:,k)),FMAdiff','*k',abs(Metrics(IndBCI,k)),FMAdiff(IndBCI)','*b',abs(Metrics(IndSham,k)),FMAdiff(IndSham)','*r');
% lsline;

% 
% 
% Electrodes16 = {'Fz','FC3','FC1','FCz','FC2','FC4','C3','C1','Cz','C2','C4','CP3','CP1','CPz','CP2','CP4'};
% i=7;
% j=8;
% figure(1004);
% plot(abs(squeeze(SMap(i,j,:))),FMAdiff','*k',abs(squeeze(SMap(i,j,IndBCI))),FMAdiff(IndBCI)','*b',abs(squeeze(SMap(i,j,IndSham))),FMAdiff(IndSham)','*r');
% lsline;
% axis([0 1.1 -2 20]);
% legend({'All','BCI','Sham'});
% ylabel('\DeltaFMA Post-Pre');
% xlabel(['Last second abs(ERSPhit-ERSPmiss) channel ' Electrodes16{i} ', Fr = ' num2str((j-1)*2+4) ' Hz']);
% [r p] = corr(abs(squeeze(SMap(i,j,:))),FMAdiff');
% text(0.4,0,['r = ' num2str(r) ', p = ' num2str(p)]);
% 
% 
% i=8;
% j=12;
% figure(1005);
% plot(abs(squeeze(SMap(i,j,:))),FMAdiff','*k',abs(squeeze(SMap(i,j,IndBCI))),FMAdiff(IndBCI)','*b',abs(squeeze(SMap(i,j,IndSham))),FMAdiff(IndSham)','*r');
% lsline;
% axis([0 1.1 -2 20]);
% legend({'All','BCI','Sham'});
% ylabel('\DeltaFMA Post-Pre');
% xlabel(['Last second abs(ERSPhit-ERSPmiss) channel ' Electrodes16{i} ', Fr = ' num2str((j-1)*2+4) ' Hz']);
% [r p] = corr(abs(squeeze(SMap(i,j,:))),FMAdiff');
% text(0.4,0,['r = ' num2str(r) ', p = ' num2str(p)]);