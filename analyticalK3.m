function [khj_33, khjme_33, kh2l_33] = analyticalK3(Vpore,Vmatrix,Vtow,km,kf,kft,beta,Pp)

ah = 1e10;
Km = km;

%% Hasselman-Johnson
lambda = sqrt(kf*kft);
lambda = kft;

A = lambda./Km - 1;
B = 1 + lambda./Km;
C = 1 - lambda./Km;
D = 1 + lambda./Km;
khj_33 = Km.*(A.*Vtow + B)./(C.*Vtow + D);

%% Hasselman-Johnson with Maxwell-Eucken correlction factor
P = Vpore;
P = (Vpore./(Vmatrix + Vpore));
%beta = 30.*Vpore;
%beta = 0;

Km = km.*((1-P)./(1+beta.*P));
% khjme_33 = kf./((1-(kf./Km)+(2*kf./ah)).*Vtow + (kf./Km));
% khjme_33 = Km.*((kf./Km-kf./ah-1).*Vtowe + (1+kf./Km+kf./ah))...
%     ./((1-kf./Km+kf./ah).*Vtowe + (1-kf./Km+kf./ah));


Vtot = Vmatrix + Vtow;
Vtow = Vtow./Vtot;

lambda = sqrt(kf*kft);
lambda = kft;

A = lambda./Km - 1;
B = 1 + lambda./Km;
C = 1 - lambda./Km;
D = 1 + lambda./Km;
khjme_33 = Km.*(A.*Vtow + B)./(C.*Vtow + D);


%khjme_33 = lambda./((1-(lambda./Km)).*Vtow + (lambda./Km));

%% H2L
P = (Vpore./(Vmatrix + Vpore));
Km = km.*((1-P)./(1+beta.*P));

A = 1 + kft./Km;
B = 1 - kft./Km;

%Pp = 0.35;
Bp = 1;

Kf = Km.*((1-Pp)./(1+Bp.*Pp)).*(1-(B./A).*Vtow)./(1+(B./A).*Vtow);

kh2l_33 = 1./(Vtow./Kf + Vmatrix./Km);

end