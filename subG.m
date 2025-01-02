function [subG_label] = subG(x,y)

if size(x)==1
    x=x';
end
if size(y)==1
    y=y';
end

Medx = median(x);
Medy = median(y);

[b,~,r,~,~] = regress(y,[ones(size(x)) x ]);

k1=b(2); 
k2=-1/b(2);

c1=Medy-k1*Medx;
c2=Medy-k2*Medx;

delta = 1;

ni=0;
while delta>=(1.5*pi)/180 && ni<=100
    %%
    for i=1:length(x)
   d01=Dis_point2line(x(i),y(i),k1,c1);
   d02=Dis_point2line(x(i),y(i),k2,c2);
   if d01<d02
       subG_label(i)=1;
   end
   if d02<d01
       subG_label(i)=2;
   end
   if d01==d02
       subG_label(i)=1.5;
   end
    end
    %%
    [b1,~,~,~,stat1] = regress(y(find(subG_label<=1.5)),...
        [ones(size(x(find(subG_label<=1.5)))) x(find(subG_label<=1.5)) ]);
    [b2,~,~,~,stat2] = regress(y(find(subG_label>=1.5)),...
        [ones(size(x(find(subG_label>=1.5)))) x(find(subG_label>=1.5)) ]);
    
    k1_new=b1(2); 
    k2_new=b2(2);

    c1_new=b1(1);
    c2_new=b2(1);

    delta=(abs(atan((k1_new-k1)/(1+k1_new*k1)))+abs(atan((k2_new-k2)/(1+k2_new*k2))))/2;
    
    k1=k1_new;
    k2=k2_new;
    c1=c1_new;
    c2=c2_new;
    
    ni=ni+1;
end