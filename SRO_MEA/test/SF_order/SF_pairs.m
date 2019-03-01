clear all;
close all;
clc;
Q=200;
FCC=textread('FCC');
ISF=textread('ISF');
load('SFE.mat');
%%%%%%%%%%%%% Monte Carlo procedure of creating random structs %%%%%%%%%%%
seed=textread('seeds');
for f=1:1:Q
    new_FCC=[];
    new_ISF=[];
    for g=1:1:length(seed(f,:));
        new_FCC(g,:)=FCC(seed(f,g),:);
        new_ISF(g,:)=ISF(seed(f,g),:);
    end
%%%%%%%%%%%%%%%% Determine SRO %%%%%%%%%%%%%%%%%%%
new_FCC(1:64,4)=1;
new_FCC(65:128,4)=2;
new_FCC(129:192,4)=3;
N=length(new_FCC);
A=10.18234;
B=10.18234;
C=24.94153;
%%%%%%%%%%%%%%%%%%%% Ghost atoms %%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_FCC(1:64,4)=1;
new_FCC(65:128,4)=2;
new_FCC(129:192,4)=3;
A=10.18234;
B=10.18234;
C=24.94153;
N_up(f)=0;
N_low(f)=0;
lower=[];
upper=[];
for i=1:1:192
   if new_FCC(i,3)>0.4*C && new_FCC(i,3)<0.43*C
       N_up(f)=N_up(f)+1;
       lower(N_up(f),1)=new_FCC(i,1);
       lower(N_up(f),2)=new_FCC(i,2);
       lower(N_up(f),3)=new_FCC(i,3);
       lower(N_up(f),4)=new_FCC(i,4);
   end
   if new_FCC(i,3)>0.49*C && new_FCC(i,3)<0.51*C
       N_low(f)=N_low(f)+1;
       upper(N_low(f),1)=new_FCC(i,1);
       upper(N_low(f),2)=new_FCC(i,2);
       upper(N_low(f),3)=new_FCC(i,3);
       upper(N_low(f),4)=new_FCC(i,4);
   end
end
o_upper=upper;
o_lower=lower;
%%%%%%%%%%%%%%%%%%%% Ghost atoms %%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(upper);
%%%%%%%%%%%%%%% Expansin along X-axis %%%%%%%%%%%%%
upper((N+1):(2*N),1)=upper(1:N,1)+A;
upper((N+1):(2*N),2)=upper(1:N,2)+0;
upper((N+1):(2*N),3)=upper(1:N,3)+0;
upper((N+1):(2*N),4)=upper(1:N,4)+0;
upper((2*N+1):(3*N),1)=upper(1:N,1)-A;
upper((2*N+1):(3*N),2)=upper(1:N,2)-0;
upper((2*N+1):(3*N),3)=upper(1:N,3)-0;
upper((2*N+1):(3*N),4)=upper(1:N,4)-0;
%%%%%%%%%%%%%%% Expansin along Y-axis %%%%%%%%%%%%%
upper((3*N+1):(6*N),1)=upper(1:(3*N),1)+B*cos(pi/3);
upper((3*N+1):(6*N),2)=upper(1:(3*N),2)+B*sin(pi/3);
upper((3*N+1):(6*N),3)=upper(1:(3*N),3)+0;
upper((3*N+1):(6*N),4)=upper(1:(3*N),4)+0;
upper((6*N+1):(9*N),1)=upper(1:(3*N),1)-B*cos(pi/3);
upper((6*N+1):(9*N),2)=upper(1:(3*N),2)-B*sin(pi/3);
upper((6*N+1):(9*N),3)=upper(1:(3*N),3)+0;
upper((6*N+1):(9*N),4)=upper(1:(3*N),4)+0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lower((N+1):(2*N),1)=lower(1:N,1)+A;
lower((N+1):(2*N),2)=lower(1:N,2)+0;
lower((N+1):(2*N),3)=lower(1:N,3)+0;
lower((N+1):(2*N),4)=lower(1:N,4)+0;
lower((2*N+1):(3*N),1)=lower(1:N,1)-A;
lower((2*N+1):(3*N),2)=lower(1:N,2)-0;
lower((2*N+1):(3*N),3)=lower(1:N,3)-0;
lower((2*N+1):(3*N),4)=lower(1:N,4)-0;
%%%%%%%%%%%%%%% Expansin along Y-axis %%%%%%%%%%%%%
lower((3*N+1):(6*N),1)=lower(1:(3*N),1)+B*cos(pi/3);
lower((3*N+1):(6*N),2)=lower(1:(3*N),2)+B*sin(pi/3);
lower((3*N+1):(6*N),3)=lower(1:(3*N),3)+0;
lower((3*N+1):(6*N),4)=lower(1:(3*N),4)+0;
lower((6*N+1):(9*N),1)=lower(1:(3*N),1)-B*cos(pi/3);
lower((6*N+1):(9*N),2)=lower(1:(3*N),2)-B*sin(pi/3);
lower((6*N+1):(9*N),3)=lower(1:(3*N),3)+0;
lower((6*N+1):(9*N),4)=lower(1:(3*N),4)+0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xu=upper(:,1);
Yu=upper(:,2);
Zu=upper(:,3);
Xl=lower(:,1);
Yl=lower(:,2);
Zl=lower(:,3);
r=2.54558;
C_a=64/192;
C_b=64/192;
C_c=64/192;
pair_type=[];
Paa(f,1)=0;
Pbb(f,1)=0;
Pcc(f,1)=0;
Pab(f,1)=0;
Pbc(f,1)=0;
Pac(f,1)=0;
for i=1:1:length(o_upper)
    k=2;
   for j=1:1:(9*length(o_upper))
       D(i,j)=sqrt((Xl(j)-Xu(i))^2+(Yl(j)-Yu(i))^2+(Zl(j)-Zu(i))^2);
       if D(i,j)>0.95*r && D(i,j)<1.05*r
           pair_type(3*i-k,1)=o_upper(i,4);
           pair_type(3*i-k,2)=lower(j,4);
           k=k-1;
       end
   end
end
for i=1:1:length(pair_type)
   if pair_type(i,1)==1 && pair_type(i,2)==1
       Paa(f,1)=Paa(f,1)+1;
   end
   if pair_type(i,1)==2 && pair_type(i,2)==2
       Pbb(f,1)=Pbb(f,1)+1;
   end
   if pair_type(i,1)==3 && pair_type(i,2)==3
       Pcc(f,1)=Pcc(f,1)+1;
   end
   if (pair_type(i,1)==1 && pair_type(i,2)==2) || (pair_type(i,1)==2 && pair_type(i,2)==1)
       Pab(f,1)=Pab(f,1)+1;
   end
   if (pair_type(i,1)==2 && pair_type(i,2)==3) || (pair_type(i,1)==3 && pair_type(i,2)==2)
       Pbc(f,1)=Pbc(f,1)+1;
   end
   if (pair_type(i,1)==1 && pair_type(i,2)==3) || (pair_type(i,1)==3 && pair_type(i,2)==1)
       Pac(f,1)=Pac(f,1)+1;
   end
end
%%%%%%%%%%%%%%%% Determine SRO %%%%%%%%%%%%%%%%%%%
new_ISF(1:64,4)=1;
new_ISF(65:128,4)=2;
new_ISF(129:192,4)=3;
N=length(new_ISF);
A=10.18234;
B=10.18234;
C=24.96646;
%%%%%%%%%%%%%%%%%%%% Ghost atoms %%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_ISF(1:64,4)=1;
new_ISF(65:128,4)=2;
new_ISF(129:192,4)=3;
A=10.18234;
B=10.18234;
C=24.96646;
N_up_isf(f)=0;
N_low_isf(f)=0;
lower_isf=[];
upper_isf=[];
for i=1:1:192
   if new_ISF(i,3)>0.4*C && new_ISF(i,3)<0.43*C
       N_up_isf(f)=N_up_isf(f)+1;
       lower_isf(N_up_isf(f),1)=new_ISF(i,1);
       lower_isf(N_up_isf(f),2)=new_ISF(i,2);
       lower_isf(N_up_isf(f),3)=new_ISF(i,3);
       lower_isf(N_up_isf(f),4)=new_ISF(i,4);
   end
   if new_ISF(i,3)>0.49*C && new_ISF(i,3)<0.51*C
       N_low_isf(f)=N_low_isf(f)+1;
       upper_isf(N_low_isf(f),1)=new_ISF(i,1);
       upper_isf(N_low_isf(f),2)=new_ISF(i,2);
       upper_isf(N_low_isf(f),3)=new_ISF(i,3);
       upper_isf(N_low_isf(f),4)=new_ISF(i,4);
   end
end
o_upper_isf=upper_isf;
o_lower_isf=lower_isf;
%%%%%%%%%%%%%%%%%%%% Ghost atoms %%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(upper_isf);
%%%%%%%%%%%%%%% Expansin along X-axis %%%%%%%%%%%%%
upper_isf((N+1):(2*N),1)=upper_isf(1:N,1)+A;
upper_isf((N+1):(2*N),2)=upper_isf(1:N,2)+0;
upper_isf((N+1):(2*N),3)=upper_isf(1:N,3)+0;
upper_isf((N+1):(2*N),4)=upper_isf(1:N,4)+0;
upper_isf((2*N+1):(3*N),1)=upper_isf(1:N,1)-A;
upper_isf((2*N+1):(3*N),2)=upper_isf(1:N,2)-0;
upper_isf((2*N+1):(3*N),3)=upper_isf(1:N,3)-0;
upper_isf((2*N+1):(3*N),4)=upper_isf(1:N,4)-0;
%%%%%%%%%%%%%%% Expansin along Y-axis %%%%%%%%%%%%%
upper_isf((3*N+1):(6*N),1)=upper_isf(1:(3*N),1)+B*cos(pi/3);
upper_isf((3*N+1):(6*N),2)=upper_isf(1:(3*N),2)+B*sin(pi/3);
upper_isf((3*N+1):(6*N),3)=upper_isf(1:(3*N),3)+0;
upper_isf((3*N+1):(6*N),4)=upper_isf(1:(3*N),4)+0;
upper_isf((6*N+1):(9*N),1)=upper_isf(1:(3*N),1)-B*cos(pi/3);
upper_isf((6*N+1):(9*N),2)=upper_isf(1:(3*N),2)-B*sin(pi/3);
upper_isf((6*N+1):(9*N),3)=upper_isf(1:(3*N),3)+0;
upper_isf((6*N+1):(9*N),4)=upper_isf(1:(3*N),4)+0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lower_isf((N+1):(2*N),1)=lower_isf(1:N,1)+A;
lower_isf((N+1):(2*N),2)=lower_isf(1:N,2)+0;
lower_isf((N+1):(2*N),3)=lower_isf(1:N,3)+0;
lower_isf((N+1):(2*N),4)=lower_isf(1:N,4)+0;
lower_isf((2*N+1):(3*N),1)=lower_isf(1:N,1)-A;
lower_isf((2*N+1):(3*N),2)=lower_isf(1:N,2)-0;
lower_isf((2*N+1):(3*N),3)=lower_isf(1:N,3)-0;
lower_isf((2*N+1):(3*N),4)=lower_isf(1:N,4)-0;
%%%%%%%%%%%%%%% Expansin along Y-axis %%%%%%%%%%%%%
lower_isf((3*N+1):(6*N),1)=lower_isf(1:(3*N),1)+B*cos(pi/3);
lower_isf((3*N+1):(6*N),2)=lower_isf(1:(3*N),2)+B*sin(pi/3);
lower_isf((3*N+1):(6*N),3)=lower_isf(1:(3*N),3)+0;
lower_isf((3*N+1):(6*N),4)=lower_isf(1:(3*N),4)+0;
lower_isf((6*N+1):(9*N),1)=lower_isf(1:(3*N),1)-B*cos(pi/3);
lower_isf((6*N+1):(9*N),2)=lower_isf(1:(3*N),2)-B*sin(pi/3);
lower_isf((6*N+1):(9*N),3)=lower_isf(1:(3*N),3)+0;
lower_isf((6*N+1):(9*N),4)=lower_isf(1:(3*N),4)+0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xu_isf=upper_isf(:,1);
Yu_isf=upper_isf(:,2);
Zu_isf=upper_isf(:,3);
Xl_isf=lower_isf(:,1);
Yl_isf=lower_isf(:,2);
Zl_isf=lower_isf(:,3);
r=2.54558;
C_a=64/192;
C_b=64/192;
C_c=64/192;
pair_type_isf=[];
Paa_isf(f,1)=0;
Pbb_isf(f,1)=0;
Pcc_isf(f,1)=0;
Pab_isf(f,1)=0;
Pbc_isf(f,1)=0;
Pac_isf(f,1)=0;
for i=1:1:length(o_upper_isf)
    k=2;
   for j=1:1:(9*length(o_upper_isf))
       D(i,j)=sqrt((Xl_isf(j)-Xu_isf(i))^2+(Yl_isf(j)-Yu_isf(i))^2+(Zl_isf(j)-Zu_isf(i))^2);
       if D(i,j)>0.95*r && D(i,j)<1.05*r
           pair_type_isf(3*i-k,1)=o_upper_isf(i,4);
           pair_type_isf(3*i-k,2)=lower_isf(j,4);
           k=k-1;
       end
   end
end
for i=1:1:length(pair_type_isf)
   if pair_type_isf(i,1)==1 && pair_type_isf(i,2)==1
   Paa_isf(f,1)=Paa_isf(f,1)+1;
   end
   if pair_type_isf(i,1)==2 && pair_type_isf(i,2)==2
       Pbb_isf(f,1)=Pbb_isf(f,1)+1;
   end
   if pair_type_isf(i,1)==3 && pair_type_isf(i,2)==3
       Pcc_isf(f,1)=Pcc_isf(f,1)+1;
   end
   if (pair_type_isf(i,1)==1 && pair_type_isf(i,2)==2) || (pair_type_isf(i,1)==2 && pair_type_isf(i,2)==1)
       Pab_isf(f,1)=Pab_isf(f,1)+1;
   end
   if (pair_type_isf(i,1)==2 && pair_type_isf(i,2)==3) || (pair_type_isf(i,1)==3 && pair_type_isf(i,2)==2)
       Pbc_isf(f,1)=Pbc_isf(f,1)+1;
   end
   if (pair_type_isf(i,1)==1 && pair_type_isf(i,2)==3) || (pair_type_isf(i,1)==3 && pair_type_isf(i,2)==1)
       Pac_isf(f,1)=Pac_isf(f,1)+1;
   end
end
end
   P_FCC(:,1)=Paa;
   P_FCC(:,2)=Pbb;
   P_FCC(:,3)=Pcc;
   P_FCC(:,4)=Pab;
   P_FCC(:,5)=Pbc;
   P_FCC(:,6)=Pac;
   P_ISF(:,1)=Paa_isf;
   P_ISF(:,2)=Pbb_isf;
   P_ISF(:,3)=Pcc_isf;
   P_ISF(:,4)=Pab_isf;
   P_ISF(:,5)=Pbc_isf;
   P_ISF(:,6)=Pac_isf;
   dpairs=P_ISF-P_FCC;