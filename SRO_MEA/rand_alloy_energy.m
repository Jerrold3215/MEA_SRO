% close all;
% clear all;
% clc;
load('analysis.mat');
data(:,1:18)=FCC;
data(:,19)=E_fcc;
data_isf(:,1:18)=ISF;
data_isf(:,19)=E_isf;
% CSRO(:,1)=data(:,2)+data(:,4)-data(:,8)-data(:,10);
CSRO(:,1)=data(:,2)+data(:,4)+data(:,8)+data(:,10);
% CSRO(:,1)=data_isf(:,2)-data_isf(:,4);
CSRO(:,2)=data(:,19);
CSRO(:,3)=SFE;
sCSRO=sortrows(CSRO);
for i=1:1:200
   A(i)=mean(sCSRO(1:i,1));
   mSFE(i)=mean(sCSRO(1:i,3));
   mean_E(i)=mean(sCSRO(1:i,2));
end
plot(A,mSFE,'-or');
hold on
% plot(A,mean_E,'-ob');
% hold on
% plot(sCSRO(:,1),mean_E,'-ob');
% hold on
% for i=1:1:200
% M1(i)=mean(FCC(1:i,1));
% M2(i)=mean(FCC(1:i,2));
% M3(i)=mean(FCC(1:i,3));
% ME_fcc(i)=mean(E_fcc(1:i,1));
% ME_isf(i)=mean(E_isf(1:i,1));
% M_SFE(i)=mean(SFE(1:i,1));
% end
% plot(M1,'-ob')
% hold on
% plot(M2,'-xr')
% hold on
% plot(M3,'-sg')
% hold on