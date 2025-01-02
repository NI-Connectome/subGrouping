function [subG_label] = subG3(x,y)

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
k2=tan(atan(k1)+pi/3);
k3=tan(atan(k1)+2*pi/3);

c1=Medy-k1*Medx;
c2=Medy-k2*Medx;
c3=Medy-k3*Medx;

delta = 1;

ni=0;
while delta>=(1.5*pi)/180 && ni<=200
    %%
    for i=1:length(x)
   d01=Dis_point2line(x(i),y(i),k1,c1);
   d02=Dis_point2line(x(i),y(i),k2,c2);
   d03=Dis_point2line(x(i),y(i),k3,c3);
   %[]=min([d01])
   if d01<=d02&&d01<=d03
       subG_label(i)=1;
   end
   if d02<=d01&&d02<=d03
       subG_label(i)=2;
   end
   if d03<=d01&&d03<=d02
       subG_label(i)=3;
   end
   if d01==d02&&d02==d03
       subG_label(i)=0;
   end
    end
    %%
    [b1,~,~,~,stat1] = regress(y(find(subG_label<=1.5)),...
        [ones(size(x(find(subG_label<=1.5)))) x(find(subG_label<=1.5)) ]);
    [b2,~,~,~,stat2] = regress(y(find(subG_label==2)),...
        [ones(size(x(find(subG_label==2)))) x(find(subG_label==2)) ]);
    [b3,~,~,~,stat3] = regress(y(find(subG_label>=2.5)),...
        [ones(size(x(find(subG_label>=2.5)))) x(find(subG_label>=2.5)) ]);
    k1_new=b1(2); 
    k2_new=b2(2);
    k3_new=b3(2);

    c1_new=b1(1);
    c2_new=b2(1);
    c3_new=b3(1);

    delta=(abs(atan((k1_new-k1)/(1+k1_new*k1)))+abs(atan((k2_new-k2)/(1+k2_new*k2)))...
        +abs(atan((k3_new-k3)/(1+k3_new*k3))))/3;
    
    k1=k1_new;
    k2=k2_new;
    k3=k3_new;
    c1=c1_new;
    c2=c2_new;
    c3=c3_new;
    ni=ni+1;
end