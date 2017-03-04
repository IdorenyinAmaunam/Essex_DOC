function printCO(Pre, Post, Group, name)

% Should add ANOVA once implemented

% Remove NaNs
IndNaN = find(isnan(Pre) | isnan(Post));
disp(['Excluded ' num2str(length(IndNaN)) ' subject due to missing data']);
Pre(IndNaN) = [];
Post(IndNaN) = [];
Group(IndNaN) = [];

pPreBCIvsPreSham = ranksum(Pre(find(strcmp(Group,'BCI'))),Pre(find(strcmp(Group,'Sham'))));
pPostBCIvsPostSham = ranksum(Post(find(strcmp(Group,'BCI'))),Post(find(strcmp(Group,'Sham'))));
pPreBCIvsPostBCI = ranksum(Pre(find(strcmp(Group,'BCI'))),Post(find(strcmp(Group,'BCI'))));
pPreShamvsPostSham = ranksum(Pre(find(strcmp(Group,'Sham'))),Post(find(strcmp(Group,'Sham'))));
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