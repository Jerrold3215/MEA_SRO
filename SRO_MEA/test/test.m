clear all;
close all;
clc;
Q=10000;
FCC=textread('FCC');
ISF=textread('ISF');
%%%%%%%%%%%%% Monte Carlo procedure of creating random structs %%%%%%%%%%%
seed=zeros(Q,192);
for f=1:1:Q
  R=randperm(192,192);
  if ismember(R,seed,'rows')==0
  seed(f,:)=R;
  end
end
seed(all(seed==0,2),:) = [];
fid=fopen('seeds','wt');
  [m,n]=size(seed);
 for i=1:1:m
    for j=1:1:n
       if j==n
         fprintf(fid,'%g\n',seed(i,j));
      else
        fprintf(fid,'%g\t',seed(i,j));
       end
    end
end
fclose(fid);
for f=1:1:Q
    new_FCC=[];
    new_ISF=[];
    for g=1:1:length(seed(f,:));
        new_FCC(g,:)=FCC(seed(f,g),:);
        new_ISF(g,:)=ISF(seed(f,g),:);
    end
%   fid=fopen(['POSCAR-',int2str(f)],'wt');
%   fprintf(fid,'%s\n','POSCAR');
%   fprintf(fid,'%s\n','1.00000000000000');
%   fprintf(fid,'%s\n','10.1823377609         0.0000000000         0.0000000000');
%   fprintf(fid,'%s\n','5.0911688805         8.8181631709         0.0000000000');
%   fprintf(fid,'%s\n','0.0000000000         0.0000000000        24.9415321350');
%   fprintf(fid,'%s\n','Co Cr Ni');
%   fprintf(fid,'%s\n','64 64 64');
%   fprintf(fid,'%s\n','Cartesian');
%   [m,n]=size(new_FCC);
%  for i=1:1:m
%     for j=1:1:n
%        if j==n
%          fprintf(fid,'%g\n',new_FCC(i,j));
%       else
%         fprintf(fid,'%g\t',new_FCC(i,j));
%        end
%     end
% end
% fclose(fid);
%   fid=fopen(['POSCAR-ISF-',int2str(f)],'wt');
%   fprintf(fid,'%s\n','POSCAR');
%   fprintf(fid,'%s\n','1.00000000000000');
%   fprintf(fid,'%s\n','10.1823377609         0.0000000000         0.0000000000');
%   fprintf(fid,'%s\n','5.0911688805         8.8181631709         0.0000000000');
%   fprintf(fid,'%s\n','0.9976617013         0.4988292166        24.9415321197');
%   fprintf(fid,'%s\n','Co Cr Ni');
%   fprintf(fid,'%s\n','64 64 64');
%   fprintf(fid,'%s\n','Cartesian');
% [m,n]=size(new_ISF);
%  for i=1:1:m
%     for j=1:1:n
%        if j==n
%          fprintf(fid,'%g\n',new_ISF(i,j));
%       else
%         fprintf(fid,'%g\t',new_ISF(i,j));
%        end
%     end
% end
% fclose(fid);
%%%%%%%%%%%%%%%% Determine SRO %%%%%%%%%%%%%%%%%%%
new_FCC(1:64,4)=1;
new_FCC(65:128,4)=2;
new_FCC(129:192,4)=3;
N=length(new_FCC);
A=10.18234;
B=10.18234;
C=24.94153;
%%%%%%%%%%%%%%%%%%%% Ghost atoms %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Expansin along X-axis %%%%%%%%%%%%%
new_FCC(1:64,4)=1;
new_FCC(65:128,4)=2;
new_FCC(129:192,4)=3;
N=length(new_FCC);
A=10.18234;
B=10.18234;
C=24.94153;
%%%%%%%%%%%%%%%%%%%% Ghost atoms %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Expansin along X-axis %%%%%%%%%%%%%
new_FCC((N+1):(2*N),1)=new_FCC(1:N,1)+A;
new_FCC((N+1):(2*N),2)=new_FCC(1:N,2)+0;
new_FCC((N+1):(2*N),3)=new_FCC(1:N,3)+0;
new_FCC((N+1):(2*N),4)=new_FCC(1:N,4)+0;
new_FCC((2*N+1):(3*N),1)=new_FCC(1:N,1)-A;
new_FCC((2*N+1):(3*N),2)=new_FCC(1:N,2)-0;
new_FCC((2*N+1):(3*N),3)=new_FCC(1:N,3)-0;
new_FCC((2*N+1):(3*N),4)=new_FCC(1:N,4)-0;
%%%%%%%%%%%%%%% Expansin along Y-axis %%%%%%%%%%%%%
new_FCC((3*N+1):(6*N),1)=new_FCC(1:(3*N),1)+B*cos(pi/3);
new_FCC((3*N+1):(6*N),2)=new_FCC(1:(3*N),2)+B*sin(pi/3);
new_FCC((3*N+1):(6*N),3)=new_FCC(1:(3*N),3)+0;
new_FCC((3*N+1):(6*N),4)=new_FCC(1:(3*N),4)+0;
new_FCC((6*N+1):(9*N),1)=new_FCC(1:(3*N),1)-B*cos(pi/3);
new_FCC((6*N+1):(9*N),2)=new_FCC(1:(3*N),2)-B*sin(pi/3);
new_FCC((6*N+1):(9*N),3)=new_FCC(1:(3*N),3)+0;
new_FCC((6*N+1):(9*N),4)=new_FCC(1:(3*N),4)+0;
%%%%%%%%%%%%%%% Expansin along Z-axis %%%%%%%%%%%%%
new_FCC((9*N+1):(18*N),1)=new_FCC(1:(9*N),1)+0;
new_FCC((9*N+1):(18*N),2)=new_FCC(1:(9*N),2)+0;
new_FCC((9*N+1):(18*N),3)=new_FCC(1:(9*N),3)+C;
new_FCC((9*N+1):(18*N),4)=new_FCC(1:(9*N),4)+0;
new_FCC((18*N+1):(27*N),1)=new_FCC(1:(9*N),1)+0;
new_FCC((18*N+1):(27*N),2)=new_FCC(1:(9*N),2)+0;
new_FCC((18*N+1):(27*N),3)=new_FCC(1:(9*N),3)-C;
new_FCC((18*N+1):(27*N),4)=new_FCC(1:(9*N),4)-0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=new_FCC(:,1);
Y=new_FCC(:,2);
Z=new_FCC(:,3);
r=2.54558;
C_a=64/192;
C_b=64/192;
C_c=64/192;
for i=1:1:N
   for j=1:1:(27*N)
       Dx(i,j)=abs(X(i)-X(j));
       Dy(i,j)=abs(Y(i)-Y(j));
       Dz(i,j)=abs(Z(i)-Z(j));
   end
end
% for i=1:1:N
%    for j=1:1:(i-1)
%        Dx(i,j)=A-abs(X(i)-X(j));
%        Dy(i,j)=B*sin(pi/3)-abs(Y(i)-Y(j));
%        Dz(i,j)=C-abs(Z(i)-Z(j));
%    end
% end
% for i=1:1:N
%    for j=1:1:N
%        if Dx(i,j)>A/1.8
%            Dx(i,j)=A-Dx(i,j);
%        end
%        if Dy(i,j)>B*sin(pi/3)/1.8
%            Dy(i,j)=B-Dy(i,j);
%        end
%        if Dz(i,j)>C/1.8
%            Dz(i,j)=C-Dz(i,j);
%        end
%    end
% end
%%%%%%%%%%%%%%%%% Count FNN_FCC SNN_FCC TNN_FCC %%%%%%%%%%%%%%%%%%%%%
for i=1:1:N
    k=1;
    m=1;
    n=1;
   for j=1:1:(27*N)
   D(i,j)=sqrt(Dx(i,j)^2+Dy(i,j)^2+Dz(i,j)^2);
   if D(i,j)> 0.95*r && D(i,j)<1.05*r
       FNN_FCC(i,k)=new_FCC(j,4);
       k=k+1;
   end
   if D(i,j)>2.8*r/2 && D(i,j)<2.9*r/2
       SNN_FCC(i,m)=new_FCC(j,4);
       m=m+1;
   end
   if D(i,j)>4.3 && D(i,j)<4.5
       TNN_FCC(i,n)=new_FCC(j,4);
       n=n+1;
   end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Determin SRO among species %%%%%%%%%%%%%
for i=1:1:64
    N_aa_FNN_FCC(i)=0;
    P_aa_FNN_FCC(i)=0;
    for j=1:1:size(FNN_FCC,2);
      if FNN_FCC(i,j)==1;
         N_aa_FNN_FCC(i)=N_aa_FNN_FCC(i)+1;
         P_aa_FNN_FCC(i)=N_aa_FNN_FCC(i)/size(FNN_FCC,2);
      end
   end
end
OP_aa_FNN_FCC=1-mean(P_aa_FNN_FCC)/C_a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=65:1:128
    N_bb_FNN_FCC(i-64)=0;
    P_bb_FNN_FCC(i-64)=0;
    for j=1:1:size(FNN_FCC,2);
      if FNN_FCC(i,j)==2;
         N_bb_FNN_FCC(i-64)=N_bb_FNN_FCC(i-64)+1;
         P_bb_FNN_FCC(i-64)=N_bb_FNN_FCC(i-64)/size(FNN_FCC,2);
      end
   end
end
OP_bb_FNN_FCC=1-mean(P_bb_FNN_FCC)/C_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_ab_FNN_FCC(i)=0;
    P_ab_FNN_FCC(i)=0;
    for j=1:1:size(FNN_FCC,2);
      if FNN_FCC(i,j)==2;
         N_ab_FNN_FCC(i)=N_ab_FNN_FCC(i)+1;
         P_ab_FNN_FCC(i)=N_ab_FNN_FCC(i)/size(FNN_FCC,2);
      end
   end
end
OP_ab_FNN_FCC=1-mean(P_ab_FNN_FCC)/C_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_ac_FNN_FCC(i)=0;
    P_ac_FNN_FCC(i)=0;
    for j=1:1:size(FNN_FCC,2);
      if FNN_FCC(i,j)==3;
         N_ac_FNN_FCC(i)=N_ac_FNN_FCC(i)+1;
         P_ac_FNN_FCC(i)=N_ac_FNN_FCC(i)/size(FNN_FCC,2);
      end
   end
end
OP_ac_FNN_FCC=1-mean(P_ac_FNN_FCC)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=129:1:192
    N_cc_FNN_FCC(i-128)=0;
    P_cc_FNN_FCC(i-128)=0;
    for j=1:1:size(FNN_FCC,2);
      if FNN_FCC(i,j)==3;
         N_cc_FNN_FCC(i-128)=N_cc_FNN_FCC(i-128)+1;
         P_cc_FNN_FCC(i-128)=N_cc_FNN_FCC(i-128)/size(FNN_FCC,2);
      end
   end
end
OP_cc_FNN_FCC=1-mean(P_cc_FNN_FCC)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=65:1:128
    N_bc_FNN_FCC(i-64)=0;
    P_bc_FNN_FCC(i-64)=0;
    for j=1:1:size(FNN_FCC,2);
      if FNN_FCC(i,j)==3;
         N_bc_FNN_FCC(i-64)=N_bc_FNN_FCC(i-64)+1;
         P_bc_FNN_FCC(i-64)=N_bc_FNN_FCC(i-64)/size(FNN_FCC,2);
      end
   end
end
OP_bc_FNN_FCC=1-mean(P_bc_FNN_FCC)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_aa_SNN_FCC(i)=0;
    P_aa_SNN_FCC(i)=0;
    for j=1:1:size(SNN_FCC,2);
      if SNN_FCC(i,j)==1;
         N_aa_SNN_FCC(i)=N_aa_SNN_FCC(i)+1;
         P_aa_SNN_FCC(i)=N_aa_SNN_FCC(i)/size(SNN_FCC,2);
      end
   end
end
OP_aa_SNN_FCC=1-mean(P_aa_SNN_FCC)/C_a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=65:1:128
    N_bb_SNN_FCC(i-64)=0;
    P_bb_SNN_FCC(i-64)=0;
    for j=1:1:size(SNN_FCC,2);
      if SNN_FCC(i,j)==2;
         N_bb_SNN_FCC(i-64)=N_bb_SNN_FCC(i-64)+1;
         P_bb_SNN_FCC(i-64)=N_bb_SNN_FCC(i-64)/size(SNN_FCC,2);
      end
   end
end
OP_bb_SNN_FCC=1-mean(P_bb_SNN_FCC)/C_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_ab_SNN_FCC(i)=0;
    P_ab_SNN_FCC(i)=0;
    for j=1:1:size(SNN_FCC,2);
      if SNN_FCC(i,j)==2;
         N_ab_SNN_FCC(i)=N_ab_SNN_FCC(i)+1;
         P_ab_SNN_FCC(i)=N_ab_SNN_FCC(i)/size(SNN_FCC,2);
      end
   end
end
OP_ab_SNN_FCC=1-mean(P_ab_SNN_FCC)/C_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_ac_SNN_FCC(i)=0;
    P_ac_SNN_FCC(i)=0;
    for j=1:1:size(SNN_FCC,2);
      if SNN_FCC(i,j)==3;
         N_ac_SNN_FCC(i)=N_ac_SNN_FCC(i)+1;
         P_ac_SNN_FCC(i)=N_ac_SNN_FCC(i)/size(SNN_FCC,2);
      end
   end
end
OP_ac_SNN_FCC=1-mean(P_ac_SNN_FCC)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=129:1:192
    N_cc_SNN_FCC(i-128)=0;
    P_cc_SNN_FCC(i-128)=0;
    for j=1:1:size(SNN_FCC,2);
      if SNN_FCC(i,j)==3;
         N_cc_SNN_FCC(i-128)=N_cc_SNN_FCC(i-128)+1;
         P_cc_SNN_FCC(i-128)=N_cc_SNN_FCC(i-128)/size(SNN_FCC,2);
      end
   end
end
OP_cc_SNN_FCC=1-mean(P_cc_SNN_FCC)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=65:1:128
    N_bc_SNN_FCC(i-64)=0;
    P_bc_SNN_FCC(i-64)=0;
    for j=1:1:size(SNN_FCC,2);
      if SNN_FCC(i,j)==3;
         N_bc_SNN_FCC(i-64)=N_bc_SNN_FCC(i-64)+1;
         P_bc_SNN_FCC(i-64)=N_bc_SNN_FCC(i-64)/size(SNN_FCC,2);
      end
   end
end
OP_bc_SNN_FCC=1-mean(P_bc_SNN_FCC)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_aa_TNN_FCC(i)=0;
    P_aa_TNN_FCC(i)=0;
    for j=1:1:size(TNN_FCC,2);
      if TNN_FCC(i,j)==1;
         N_aa_TNN_FCC(i)=N_aa_TNN_FCC(i)+1;
         P_aa_TNN_FCC(i)=N_aa_TNN_FCC(i)/size(TNN_FCC,2);
      end
   end
end
OP_aa_TNN_FCC=1-mean(P_aa_TNN_FCC)/C_a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=65:1:128
    N_bb_TNN_FCC(i-64)=0;
    P_bb_TNN_FCC(i-64)=0;
    for j=1:1:size(TNN_FCC,2);
      if TNN_FCC(i,j)==2;
         N_bb_TNN_FCC(i-64)=N_bb_TNN_FCC(i-64)+1;
         P_bb_TNN_FCC(i-64)=N_bb_TNN_FCC(i-64)/size(TNN_FCC,2);
      end
   end
end
OP_bb_TNN_FCC=1-mean(P_bb_TNN_FCC)/C_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_ab_TNN_FCC(i)=0;
    P_ab_TNN_FCC(i)=0;
    for j=1:1:size(TNN_FCC,2);
      if TNN_FCC(i,j)==2;
         N_ab_TNN_FCC(i)=N_ab_TNN_FCC(i)+1;
         P_ab_TNN_FCC(i)=N_ab_TNN_FCC(i)/size(TNN_FCC,2);
      end
   end
end
OP_ab_TNN_FCC=1-mean(P_ab_TNN_FCC)/C_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_ac_TNN_FCC(i)=0;
    P_ac_TNN_FCC(i)=0;
    for j=1:1:size(TNN_FCC,2);
      if TNN_FCC(i,j)==3;
         N_ac_TNN_FCC(i)=N_ac_TNN_FCC(i)+1;
         P_ac_TNN_FCC(i)=N_ac_TNN_FCC(i)/size(TNN_FCC,2);
      end
   end
end
OP_ac_TNN_FCC=1-mean(P_ac_TNN_FCC)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=129:1:192
    N_cc_TNN_FCC(i-128)=0;
    P_cc_TNN_FCC(i-128)=0;
    for j=1:1:size(TNN_FCC,2);
      if TNN_FCC(i,j)==3;
         N_cc_TNN_FCC(i-128)=N_cc_TNN_FCC(i-128)+1;
         P_cc_TNN_FCC(i-128)=N_cc_TNN_FCC(i-128)/size(TNN_FCC,2);
      end
   end
end
OP_cc_TNN_FCC=1-mean(P_cc_TNN_FCC)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=65:1:128
    N_bc_TNN_FCC(i-64)=0;
    P_bc_TNN_FCC(i-64)=0;
    for j=1:1:size(TNN_FCC,2);
      if TNN_FCC(i,j)==3;
         N_bc_TNN_FCC(i-64)=N_bc_TNN_FCC(i-64)+1;
         P_bc_TNN_FCC(i-64)=N_bc_TNN_FCC(i-64)/size(TNN_FCC,2);
      end
   end
end
OP_bc_TNN_FCC=1-mean(P_bc_TNN_FCC)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_ISF(1:64,4)=1;
new_ISF(65:128,4)=2;
new_ISF(129:192,4)=3;
N=length(new_ISF);
A=10.18234;
B=10.18234;
C=24.94153;
%%%%%%%%%%%%%%%%%%%% Ghost atoms %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Expansin along X-axis %%%%%%%%%%%%%
new_ISF((N+1):(2*N),1)=new_ISF(1:N,1)+A;
new_ISF((N+1):(2*N),2)=new_ISF(1:N,2)+0;
new_ISF((N+1):(2*N),3)=new_ISF(1:N,3)+0;
new_ISF((N+1):(2*N),4)=new_ISF(1:N,4)+0;
new_ISF((2*N+1):(3*N),1)=new_ISF(1:N,1)-A;
new_ISF((2*N+1):(3*N),2)=new_ISF(1:N,2)-0;
new_ISF((2*N+1):(3*N),3)=new_ISF(1:N,3)-0;
new_ISF((2*N+1):(3*N),4)=new_ISF(1:N,4)-0;
%%%%%%%%%%%%%%% Expansin along Y-axis %%%%%%%%%%%%%
new_ISF((3*N+1):(6*N),1)=new_ISF(1:(3*N),1)+B*cos(pi/3);
new_ISF((3*N+1):(6*N),2)=new_ISF(1:(3*N),2)+B*sin(pi/3);
new_ISF((3*N+1):(6*N),3)=new_ISF(1:(3*N),3)+0;
new_ISF((3*N+1):(6*N),4)=new_ISF(1:(3*N),4)+0;
new_ISF((6*N+1):(9*N),1)=new_ISF(1:(3*N),1)-B*cos(pi/3);
new_ISF((6*N+1):(9*N),2)=new_ISF(1:(3*N),2)-B*sin(pi/3);
new_ISF((6*N+1):(9*N),3)=new_ISF(1:(3*N),3)+0;
new_ISF((6*N+1):(9*N),4)=new_ISF(1:(3*N),4)+0;
%%%%%%%%%%%%%%% Expansin along Z-axis %%%%%%%%%%%%%
new_ISF((9*N+1):(18*N),1)=new_ISF(1:(9*N),1)+0;
new_ISF((9*N+1):(18*N),2)=new_ISF(1:(9*N),2)+0;
new_ISF((9*N+1):(18*N),3)=new_ISF(1:(9*N),3)+C;
new_ISF((9*N+1):(18*N),4)=new_ISF(1:(9*N),4)+0;
new_ISF((18*N+1):(27*N),1)=new_ISF(1:(9*N),1)+0;
new_ISF((18*N+1):(27*N),2)=new_ISF(1:(9*N),2)+0;
new_ISF((18*N+1):(27*N),3)=new_ISF(1:(9*N),3)-C;
new_ISF((18*N+1):(27*N),4)=new_ISF(1:(9*N),4)-0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=new_ISF(:,1);
Y=new_ISF(:,2);
Z=new_ISF(:,3);
r=2.54558;
C_a=64/192;
C_b=64/192;
C_c=64/192;
for i=1:1:N
   for j=1:1:(27*N)
       Dx(i,j)=abs(X(i)-X(j));
       Dy(i,j)=abs(Y(i)-Y(j));
       Dz(i,j)=abs(Z(i)-Z(j));
   end
end
% for i=1:1:N
%    for j=1:1:(i-1)
%        Dx(i,j)=A-abs(X(i)-X(j));
%        Dy(i,j)=B*sin(pi/3)-abs(Y(i)-Y(j));
%        Dz(i,j)=C-abs(Z(i)-Z(j));
%    end
% end
% for i=1:1:N
%    for j=1:1:N
%        if Dx(i,j)>A/1.8
%            Dx(i,j)=A-Dx(i,j);
%        end
%        if Dy(i,j)>B*sin(pi/3)/1.8
%            Dy(i,j)=B-Dy(i,j);
%        end
%        if Dz(i,j)>C/1.8
%            Dz(i,j)=C-Dz(i,j);
%        end
%    end
% end
%%%%%%%%%%%%%%%%% Count FNN_ISF SNN_ISF TNN_ISF %%%%%%%%%%%%%%%%%%%%%
for i=1:1:N
    k=1;
    m=1;
    n=1;
   for j=1:1:(27*N)
   D(i,j)=sqrt(Dx(i,j)^2+Dy(i,j)^2+Dz(i,j)^2);
   if D(i,j)> 0.95*r && D(i,j)<1.05*r
       FNN_ISF(i,k)=new_ISF(j,4);
       k=k+1;
   end
   if D(i,j)>2.8*r/2 && D(i,j)<2.9*r/2
       SNN_ISF(i,m)=new_ISF(j,4);
       m=m+1;
   end
   if D(i,j)>4.3 && D(i,j)<4.5
       TNN_ISF(i,n)=new_ISF(j,4);
       n=n+1;
   end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Determin SRO among species %%%%%%%%%%%%%
for i=1:1:64
    N_aa_FNN_ISF(i)=0;
    P_aa_FNN_ISF(i)=0;
    for j=1:1:size(FNN_ISF,2);
      if FNN_ISF(i,j)==1;
         N_aa_FNN_ISF(i)=N_aa_FNN_ISF(i)+1;
         P_aa_FNN_ISF(i)=N_aa_FNN_ISF(i)/size(FNN_ISF,2);
      end
   end
end
OP_aa_FNN_ISF=1-mean(P_aa_FNN_ISF)/C_a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=65:1:128
    N_bb_FNN_ISF(i-64)=0;
    P_bb_FNN_ISF(i-64)=0;
    for j=1:1:size(FNN_ISF,2);
      if FNN_ISF(i,j)==2;
         N_bb_FNN_ISF(i-64)=N_bb_FNN_ISF(i-64)+1;
         P_bb_FNN_ISF(i-64)=N_bb_FNN_ISF(i-64)/size(FNN_ISF,2);
      end
   end
end
OP_bb_FNN_ISF=1-mean(P_bb_FNN_ISF)/C_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_ab_FNN_ISF(i)=0;
    P_ab_FNN_ISF(i)=0;
    for j=1:1:size(FNN_ISF,2);
      if FNN_ISF(i,j)==2;
         N_ab_FNN_ISF(i)=N_ab_FNN_ISF(i)+1;
         P_ab_FNN_ISF(i)=N_ab_FNN_ISF(i)/size(FNN_ISF,2);
      end
   end
end
OP_ab_FNN_ISF=1-mean(P_ab_FNN_ISF)/C_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_ac_FNN_ISF(i)=0;
    P_ac_FNN_ISF(i)=0;
    for j=1:1:size(FNN_ISF,2);
      if FNN_ISF(i,j)==3;
         N_ac_FNN_ISF(i)=N_ac_FNN_ISF(i)+1;
         P_ac_FNN_ISF(i)=N_ac_FNN_ISF(i)/size(FNN_ISF,2);
      end
   end
end
OP_ac_FNN_ISF=1-mean(P_ac_FNN_ISF)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=129:1:192
    N_cc_FNN_ISF(i-128)=0;
    P_cc_FNN_ISF(i-128)=0;
    for j=1:1:size(FNN_ISF,2);
      if FNN_ISF(i,j)==3;
         N_cc_FNN_ISF(i-128)=N_cc_FNN_ISF(i-128)+1;
         P_cc_FNN_ISF(i-128)=N_cc_FNN_ISF(i-128)/size(FNN_ISF,2);
      end
   end
end
OP_cc_FNN_ISF=1-mean(P_cc_FNN_ISF)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=65:1:128
    N_bc_FNN_ISF(i-64)=0;
    P_bc_FNN_ISF(i-64)=0;
    for j=1:1:size(FNN_ISF,2);
      if FNN_ISF(i,j)==3;
         N_bc_FNN_ISF(i-64)=N_bc_FNN_ISF(i-64)+1;
         P_bc_FNN_ISF(i-64)=N_bc_FNN_ISF(i-64)/size(FNN_ISF,2);
      end
   end
end
OP_bc_FNN_ISF=1-mean(P_bc_FNN_ISF)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_aa_SNN_ISF(i)=0;
    P_aa_SNN_ISF(i)=0;
    for j=1:1:size(SNN_ISF,2);
      if SNN_ISF(i,j)==1;
         N_aa_SNN_ISF(i)=N_aa_SNN_ISF(i)+1;
         P_aa_SNN_ISF(i)=N_aa_SNN_ISF(i)/size(SNN_ISF,2);
      end
   end
end
OP_aa_SNN_ISF=1-mean(P_aa_SNN_ISF)/C_a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=65:1:128
    N_bb_SNN_ISF(i-64)=0;
    P_bb_SNN_ISF(i-64)=0;
    for j=1:1:size(SNN_ISF,2);
      if SNN_ISF(i,j)==2;
         N_bb_SNN_ISF(i-64)=N_bb_SNN_ISF(i-64)+1;
         P_bb_SNN_ISF(i-64)=N_bb_SNN_ISF(i-64)/size(SNN_ISF,2);
      end
   end
end
OP_bb_SNN_ISF=1-mean(P_bb_SNN_ISF)/C_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_ab_SNN_ISF(i)=0;
    P_ab_SNN_ISF(i)=0;
    for j=1:1:size(SNN_ISF,2);
      if SNN_ISF(i,j)==2;
         N_ab_SNN_ISF(i)=N_ab_SNN_ISF(i)+1;
         P_ab_SNN_ISF(i)=N_ab_SNN_ISF(i)/size(SNN_ISF,2);
      end
   end
end
OP_ab_SNN_ISF=1-mean(P_ab_SNN_ISF)/C_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_ac_SNN_ISF(i)=0;
    P_ac_SNN_ISF(i)=0;
    for j=1:1:size(SNN_ISF,2);
      if SNN_ISF(i,j)==3;
         N_ac_SNN_ISF(i)=N_ac_SNN_ISF(i)+1;
         P_ac_SNN_ISF(i)=N_ac_SNN_ISF(i)/size(SNN_ISF,2);
      end
   end
end
OP_ac_SNN_ISF=1-mean(P_ac_SNN_ISF)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=129:1:192
    N_cc_SNN_ISF(i-128)=0;
    P_cc_SNN_ISF(i-128)=0;
    for j=1:1:size(SNN_ISF,2);
      if SNN_ISF(i,j)==3;
         N_cc_SNN_ISF(i-128)=N_cc_SNN_ISF(i-128)+1;
         P_cc_SNN_ISF(i-128)=N_cc_SNN_ISF(i-128)/size(SNN_ISF,2);
      end
   end
end
OP_cc_SNN_ISF=1-mean(P_cc_SNN_ISF)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=65:1:128
    N_bc_SNN_ISF(i-64)=0;
    P_bc_SNN_ISF(i-64)=0;
    for j=1:1:size(SNN_ISF,2);
      if SNN_ISF(i,j)==3;
         N_bc_SNN_ISF(i-64)=N_bc_SNN_ISF(i-64)+1;
         P_bc_SNN_ISF(i-64)=N_bc_SNN_ISF(i-64)/size(SNN_ISF,2);
      end
   end
end
OP_bc_SNN_ISF=1-mean(P_bc_SNN_ISF)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_aa_TNN_ISF(i)=0;
    P_aa_TNN_ISF(i)=0;
    for j=1:1:size(TNN_ISF,2);
      if TNN_ISF(i,j)==1;
         N_aa_TNN_ISF(i)=N_aa_TNN_ISF(i)+1;
         P_aa_TNN_ISF(i)=N_aa_TNN_ISF(i)/size(TNN_ISF,2);
      end
   end
end
OP_aa_TNN_ISF=1-mean(P_aa_TNN_ISF)/C_a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=65:1:128
    N_bb_TNN_ISF(i-64)=0;
    P_bb_TNN_ISF(i-64)=0;
    for j=1:1:size(TNN_ISF,2);
      if TNN_ISF(i,j)==2;
         N_bb_TNN_ISF(i-64)=N_bb_TNN_ISF(i-64)+1;
         P_bb_TNN_ISF(i-64)=N_bb_TNN_ISF(i-64)/size(TNN_ISF,2);
      end
   end
end
OP_bb_TNN_ISF=1-mean(P_bb_TNN_ISF)/C_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_ab_TNN_ISF(i)=0;
    P_ab_TNN_ISF(i)=0;
    for j=1:1:size(TNN_ISF,2);
      if TNN_ISF(i,j)==2;
         N_ab_TNN_ISF(i)=N_ab_TNN_ISF(i)+1;
         P_ab_TNN_ISF(i)=N_ab_TNN_ISF(i)/size(TNN_ISF,2);
      end
   end
end
OP_ab_TNN_ISF=1-mean(P_ab_TNN_ISF)/C_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:64
    N_ac_TNN_ISF(i)=0;
    P_ac_TNN_ISF(i)=0;
    for j=1:1:size(TNN_ISF,2);
      if TNN_ISF(i,j)==3;
         N_ac_TNN_ISF(i)=N_ac_TNN_ISF(i)+1;
         P_ac_TNN_ISF(i)=N_ac_TNN_ISF(i)/size(TNN_ISF,2);
      end
   end
end
OP_ac_TNN_ISF=1-mean(P_ac_TNN_ISF)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=129:1:192
    N_cc_TNN_ISF(i-128)=0;
    P_cc_TNN_ISF(i-128)=0;
    for j=1:1:size(TNN_ISF,2);
      if TNN_ISF(i,j)==3;
         N_cc_TNN_ISF(i-128)=N_cc_TNN_ISF(i-128)+1;
         P_cc_TNN_ISF(i-128)=N_cc_TNN_ISF(i-128)/size(TNN_ISF,2);
      end
   end
end
OP_cc_TNN_ISF=1-mean(P_cc_TNN_ISF)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=65:1:128
    N_bc_TNN_ISF(i-64)=0;
    P_bc_TNN_ISF(i-64)=0;
    for j=1:1:size(TNN_ISF,2);
      if TNN_ISF(i,j)==3;
         N_bc_TNN_ISF(i-64)=N_bc_TNN_ISF(i-64)+1;
         P_bc_TNN_ISF(i-64)=N_bc_TNN_ISF(i-64)/size(TNN_ISF,2);
      end
   end
end
OP_bc_TNN_ISF=1-mean(P_bc_TNN_ISF)/C_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Record Order Parameters of each struct %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Sequence is 'aa bb cc ab bc ac' %%%%%%%%%%%%%%%%%%%%%%
OP_FCC(f,1)=OP_aa_FNN_FCC;
OP_FCC(f,2)=OP_bb_FNN_FCC;
OP_FCC(f,3)=OP_cc_FNN_FCC;
OP_FCC(f,4)=OP_ab_FNN_FCC;
OP_FCC(f,5)=OP_bc_FNN_FCC;
OP_FCC(f,6)=OP_ac_FNN_FCC;
OP_FCC(f,7)=OP_aa_SNN_FCC;
OP_FCC(f,8)=OP_bb_SNN_FCC;
OP_FCC(f,9)=OP_cc_SNN_FCC;
OP_FCC(f,10)=OP_ab_SNN_FCC;
OP_FCC(f,11)=OP_bc_SNN_FCC;
OP_FCC(f,12)=OP_ac_SNN_FCC;
OP_FCC(f,13)=OP_aa_TNN_FCC;
OP_FCC(f,14)=OP_bb_TNN_FCC;
OP_FCC(f,15)=OP_cc_TNN_FCC;
OP_FCC(f,16)=OP_ab_TNN_FCC;
OP_FCC(f,17)=OP_bc_TNN_FCC;
OP_FCC(f,18)=OP_ac_TNN_FCC;
OP_ISF(f,1)=OP_aa_FNN_ISF;
OP_ISF(f,2)=OP_bb_FNN_ISF;
OP_ISF(f,3)=OP_cc_FNN_ISF;
OP_ISF(f,4)=OP_ab_FNN_ISF;
OP_ISF(f,5)=OP_bc_FNN_ISF;
OP_ISF(f,6)=OP_ac_FNN_ISF;
OP_ISF(f,7)=OP_aa_SNN_ISF;
OP_ISF(f,8)=OP_bb_SNN_ISF;
OP_ISF(f,9)=OP_cc_SNN_ISF;
OP_ISF(f,10)=OP_ab_SNN_ISF;
OP_ISF(f,11)=OP_bc_SNN_ISF;
OP_ISF(f,12)=OP_ac_SNN_ISF;
OP_ISF(f,13)=OP_aa_TNN_ISF;
OP_ISF(f,14)=OP_bb_TNN_ISF;
OP_ISF(f,15)=OP_cc_TNN_ISF;
OP_ISF(f,16)=OP_ab_TNN_ISF;
OP_ISF(f,17)=OP_bc_TNN_ISF;
OP_ISF(f,18)=OP_ac_TNN_ISF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fid=fopen('FCC_SRO','wt');
  fprintf(fid,'%s\n','Faa Fbb Fcc Fab Fbc Fac Saa Sbb Scc Sab Sbc Sac Taa Tbb Tcc Tab Tbc Tac');
  [m,n]=size(OP_FCC);
 for i=1:1:m
    for j=1:1:n
       if j==n
         fprintf(fid,'%g\n',OP_FCC(i,j));
      else
        fprintf(fid,'%g\t',OP_FCC(i,j));
       end
    end
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fid=fopen('ISF_SRO','wt');
  fprintf(fid,'%s\n','Faa Fbb Fcc Fab Fbc Fac Saa Sbb Scc Sab Sbc Sac Taa Tbb Tcc Tab Tbc Tac');
  [m,n]=size(OP_ISF);
 for i=1:1:m
    for j=1:1:n
       if j==n
         fprintf(fid,'%g\n',OP_ISF(i,j));
      else
        fprintf(fid,'%g\t',OP_ISF(i,j));
       end
    end
end
fclose(fid);
end
for i=1:1:Q
M_fcc_1(i)=mean(OP_FCC(1:i,1));
M_fcc_2(i)=mean(OP_FCC(1:i,2));
M_fcc_3(i)=mean(OP_FCC(1:i,3));
M_isf_1(i)=mean(OP_ISF(1:i,1));
M_isf_2(i)=mean(OP_ISF(1:i,2));
M_isf_3(i)=mean(OP_ISF(1:i,3));
end
plot(M_fcc_1,'-ob')
hold on
plot(M_fcc_2,'-xr')
hold on
plot(M_fcc_3,'-sg')
hold on
% plot(ME_fcc,'-sr')
% hold on
% plot(ME_isf,'-bd')
% hold on