function [khj_33, khjme_33] = analyticalK1(Vpore,Vmatrix,Vtow,km,kf)

beta = 6.5;
P = Vpore;
ah = (mean([0.985e-3,0.115e-3])).*(km./(0.2e-6));

%% Hasselman-Johnson
Km = km;
khj_33 = Km.*((kf./Km-kf./ah-1).*Vtow + (1+kf./Km+kf./ah))...
    ./((1-kf./Km+kf./ah).*Vtow + (1-kf./Km+kf./ah));

%% Hasselman-Johnson with Maxwell-Eucken correlction factor
Km = km.*((1-P)./(1+beta.*P));
khjme_33 = Km.*((kf./Km-kf./ah-1).*Vtow + (1+kf./Km+kf./ah))...
    ./((1-kf./Km+kf./ah).*Vtow + (1-kf./Km+kf./ah));

end