function nERSPMap =  loweResERSP(ERSPMap, LocReg, FrReg)

RegIpsi = LocReg{1};
RegMed = LocReg{2};
RegContra = LocReg{3};

FrMu = FrReg{1};
FrLowBeta = FrReg{2};
FrHighBeta = FrReg{3};


nERSPMap(:,1,1) = nanmean((nanmean(ERSPMap(:,RegIpsi,FrMu),2)),3);
nERSPMap(:,1,2) = nanmean((nanmean(ERSPMap(:,RegIpsi,FrLowBeta),2)),3);
nERSPMap(:,1,3) = nanmean((nanmean(ERSPMap(:,RegIpsi,FrHighBeta),2)),3);
nERSPMap(:,2,1) = nanmean((nanmean(ERSPMap(:,RegMed,FrMu),2)),3);
nERSPMap(:,2,2) = nanmean((nanmean(ERSPMap(:,RegMed,FrLowBeta),2)),3);
nERSPMap(:,2,3) = nanmean((nanmean(ERSPMap(:,RegMed,FrHighBeta),2)),3);
nERSPMap(:,3,1) = nanmean((nanmean(ERSPMap(:,RegContra,FrMu),2)),3);
nERSPMap(:,3,2) = nanmean((nanmean(ERSPMap(:,RegContra,FrLowBeta),2)),3);
nERSPMap(:,3,3) = nanmean((nanmean(ERSPMap(:,RegContra,FrHighBeta),2)),3);
