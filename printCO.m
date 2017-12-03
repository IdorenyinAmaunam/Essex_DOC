function printCO(Pre, Post, Group, name)

Color = {[255 40 40]/255,[102 184 184]/255};
% Should add ANOVA once implemented

% Remove NaNs
IndNaN = find(isnan(Pre) | isnan(Post));
disp(['Excluded ' num2str(length(IndNaN)) ' subject due to missing data']);
Pre(IndNaN) = [];
Post(IndNaN) = [];
Group(IndNaN) = [];

% pPreBCIvsPreSham = ranksum(Pre(find(strcmp(Group,'BCI'))),Pre(find(strcmp(Group,'Sham'))));
% pPostBCIvsPostSham = ranksum(Post(find(strcmp(Group,'BCI'))),Post(find(strcmp(Group,'Sham'))));
% pPreBCIvsPostBCI = ranksum(Pre(find(strcmp(Group,'BCI'))),Post(find(strcmp(Group,'BCI'))));
% pPreShamvsPostSham = ranksum(Pre(find(strcmp(Group,'Sham'))),Post(find(strcmp(Group,'Sham'))));

[~, pPreBCIvsPreSham] = ttest2(Pre(find(strcmp(Group,'BCI'))),Pre(find(strcmp(Group,'Sham'))));
[~, pPostBCIvsPostSham] = ttest2(Post(find(strcmp(Group,'BCI'))),Post(find(strcmp(Group,'Sham'))));
[~, pPreBCIvsPostBCI] = ttest(Pre(find(strcmp(Group,'BCI'))),Post(find(strcmp(Group,'BCI'))));
[~, pPreShamvsPostSham] = ttest(Pre(find(strcmp(Group,'Sham'))),Post(find(strcmp(Group,'Sham'))));

MPreBCI = mean(Pre(find(strcmp(Group,'BCI'))));
MPostBCI = mean(Post(find(strcmp(Group,'BCI'))));
MPreSham = mean(Pre(find(strcmp(Group,'Sham'))));
MPostSham = mean(Post(find(strcmp(Group,'Sham'))));
SPreBCI = std2(Pre(find(strcmp(Group,'BCI'))));
SPostBCI = std2(Post(find(strcmp(Group,'BCI'))));
SPreSham = std2(Pre(find(strcmp(Group,'Sham'))));
SPostSham = std2(Post(find(strcmp(Group,'Sham'))));

Diff = Post-Pre;
pDiffBCIvsDiffSham = ranksum(Diff(find(strcmp(Group,'BCI'))),Diff(find(strcmp(Group,'Sham'))));
MDiffBCI = mean(Diff(find(strcmp(Group,'BCI'))));
MDiffSham = mean(Diff(find(strcmp(Group,'Sham'))));
SDiffBCI = std2(Diff(find(strcmp(Group,'BCI'))));
SDiffSham = std2(Diff(find(strcmp(Group,'Sham'))));

% Print all results
disp([name ' Pre/BCI vs ' name ' Pre/Sham : ' num2str(MPreBCI) ' +/- ' num2str(SPreBCI) '  vs  ' num2str(MPreSham) ' +/- ' num2str(SPreSham) ', p = ' num2str(pPreBCIvsPreSham)]);
disp([name ' Post/BCI vs ' name ' Post/Sham : ' num2str(MPostBCI) ' +/- ' num2str(SPostBCI) '  vs  ' num2str(MPostSham) ' +/- ' num2str(SPostSham) ', p = ' num2str(pPostBCIvsPostSham)]);
disp([name ' Pre/BCI vs ' name ' Post/BCI : ' num2str(MPreBCI) ' +/- ' num2str(SPreBCI) '  vs  ' num2str(MPostBCI) ' +/- ' num2str(SPostBCI) ', p = ' num2str(pPreBCIvsPostBCI)]);
disp([name ' Pre/Sham vs ' name ' Post/Sham : ' num2str(MPreSham) ' +/- ' num2str(SPreSham) '  vs  ' num2str(MPostSham) ' +/- ' num2str(SPostSham) ', p = ' num2str(pPreShamvsPostSham)]);
disp([name ' Difference BCI vs ' name ' Difference Sham : ' num2str(MDiffBCI) ' +/- ' num2str(SDiffBCI) '  vs  ' num2str(MDiffSham) ' +/- ' num2str(SDiffSham) ', p = ' num2str(pDiffBCIvsDiffSham)]);

% Plot all results
figure();
bar(1,MPreBCI,'FaceColor',Color{1});
hold on;
bar(2,MPreSham,'FaceColor',Color{2});
bar(4,MPostBCI,'FaceColor',Color{1});
bar(5,MPostSham,'FaceColor',Color{2});
bar(8,MDiffBCI,'FaceColor',Color{1});
bar(9,MDiffSham,'FaceColor',Color{2});
errorbar(1,MPreBCI,0,SPreBCI,'k','LineWidth',3);
errorbar(2,MPreSham,0,SPreSham,'k','LineWidth',3);
errorbar(4,MPostBCI,0,SPostBCI,'k','LineWidth',3);
errorbar(5,MPostSham,0,SPostSham,'k','LineWidth',3);
errorbar(8,MDiffBCI,0,SDiffBCI,'k','LineWidth',3);
errorbar(9,MDiffSham,0,SDiffSham,'k','LineWidth',3);
legend({'BCI','Sham'});
ylabel(name);
set(gca, 'XTick', [1.5 4.5 8.5]);
set(gca, 'XTickLabel', {'Pre-intervention','Post-intervention','Difference'},'FontSize',20);

% Build sigstar strings
strPairs = '{';
strVals = '[';
if(pPreBCIvsPreSham < 0.05)
    strPairs = [strPairs '[1,2],'];
    strVals = [strVals 'pPreBCIvsPreSham '];
end
if(pPostBCIvsPostSham < 0.05)
    strPairs = [strPairs '[4,5],'];
    strVals = [strVals 'pPostBCIvsPostSham '];    
end
if(pPreBCIvsPostBCI < 0.05)
    strPairs = [strPairs '[1,4],'];
    strVals = [strVals 'pPreBCIvsPostBCI '];    
end
if(pPreShamvsPostSham < 0.05)
    strPairs = [strPairs '[2,5],'];
    strVals = [strVals 'pPreShamvsPostSham '];    
end
if(pDiffBCIvsDiffSham < 0.05)
    strPairs = [strPairs '[8,9],'];
    strVals = [strVals 'pDiffBCIvsDiffSham '];    
end
strPairs = [strPairs '}'];
strVals = [strVals ']'];

strAll = ['sigstar(' strPairs ',' strVals ');'];
eval(strAll)
title(name);
hold off;
drawnow;