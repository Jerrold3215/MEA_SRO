% for i=1:1:100
%     N(i)=rand;
% end
% for i=1:1:100
% plot(i,mean(N(1:i)),'or')
% hold on
% end


for i=1:1:200
M1(i)=mean(FCC(1:i,1));
M2(i)=mean(FCC(1:i,2));
M3(i)=mean(FCC(1:i,3));
ME_fcc(i)=mean(E_fcc(1:i,1));
ME_isf(i)=mean(E_isf(1:i,1));
M_SFE(i)=mean(SFE(1:i,1));
end
% plot(M1,'-ob')
% hold on
% plot(M2,'-xr')
% hold on
% plot(M3,'-sg')
% hold on
% plot(ME_fcc,'-sr')
% hold on
% plot(ME_isf,'-bd')
% hold on
plot(M_SFE,'-or')
hold on