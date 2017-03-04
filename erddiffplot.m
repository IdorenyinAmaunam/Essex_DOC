function erddiffplot(PreERSP,PostERSP,IndBCI,IndSham,SelFeatCh,SelFeatFr)

% Select feature to show in bars
figure;

bar(1,nanmean(PreERSP(IndBCI,SelFeatCh,SelFeatFr),1),'b');hold on;
bar(2,nanmean(PreERSP(IndSham,SelFeatCh,SelFeatFr),1),'r');
bar(4,nanmean(PostERSP(IndBCI,SelFeatCh,SelFeatFr),1),'b');
bar(5,nanmean(PostERSP(IndSham,SelFeatCh,SelFeatFr),1),'r');
errorbar(1,nanmean(PreERSP(IndBCI,SelFeatCh,SelFeatFr),1),nanstd(PreERSP(IndBCI,SelFeatCh,SelFeatFr),1),'b','LineWidth',3);
errorbar(2,nanmean(PreERSP(IndSham,SelFeatCh,SelFeatFr),1),nanstd(PreERSP(IndSham,SelFeatCh,SelFeatFr),1),'r','LineWidth',3);
errorbar(4,nanmean(PostERSP(IndBCI,SelFeatCh,SelFeatFr),1),nanstd(PostERSP(IndBCI,SelFeatCh,SelFeatFr),1),'b','LineWidth',3);
errorbar(5,nanmean(PostERSP(IndSham,SelFeatCh,SelFeatFr),1),nanstd(PostERSP(IndSham,SelFeatCh,SelFeatFr),1),'r','LineWidth',3);hold off;
ylabel('ERSP (z-score)');
set(gca,'XTick',1.5:3:5);
set(gca,'XTickLabel',{'Pre-intervention','Post-intervention'});
legend({'BCI','Sham'});
set(gca,'LineWidth',3);
set(gca,'FontSize',18);
