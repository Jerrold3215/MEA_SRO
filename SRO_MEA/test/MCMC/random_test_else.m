for i=1:1:100
   r(i)=rand
   if r(i)<0.3
      B(i)=0
   elseif r(i)>0.7
       B(i)=1
   elseif r(i)>=0.3 && r(i)<=0.7
       B(i)=999
   end
end