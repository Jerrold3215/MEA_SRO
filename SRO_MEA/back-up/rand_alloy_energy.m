% for i=1:1:100
%     N(i)=rand;
% end
% for i=1:1:100
% plot(i,mean(N(1:i)),'or')
% hold on
% end


for i=1:1:200
plot(i,mean(FCC(1:i,2)),'xr')
plot(i,mean(ISF(1:i,2)),'xg')
% plot(i,mean(FCC(1:i,6)),'xb')
hold on
end