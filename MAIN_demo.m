%load demo.mat

clc;clear;
Nsubj=200;          % the number of selected subjects
IndHomo=1:200;  % index of subjects if select a subset of the datasets,with homo information; 
Nr=100;             % the number of permutation test for subgrouping;
Np=20;              % the number of features, for a 90x90 matrix, is 4005
Ntop=10;           % select the top features

t1=randi(50,100,1);
y1=t1+20*randn(size(t1))-25;
t2=randi(50,100,1);
y2=-t2+20*randn(size(t2))+25;
X=repmat([t1;t2],[1,20]);
X=X+randn(size(X));

Y=[y1;y2];
%% connection/feature level subgrouping 
    for m=1:Nr+1;
        % the first is real data, others are permutation tests
        for n=1:Np%
        %parfor n=1:4005;
            % running across all featues 
        x0=X(IndHomo,n);
        % [~,~,y,~,~] = regress(Y(IndHomo,1),[ones([length(IndHomo) 1]) covs]); 
        [~,~,y,~,~] = regress(Y(IndHomo,1),[ones([length(IndHomo) 1]) ]); 
           if m>1
           IndP = randperm(Nsubj);
           y = y(IndP);
           end
         %%
            %[~,~,x,~,~] = regress(x0,[ones([length(IndHomo) 1]) covs]);
            [~,~,x,~,~] = regress(x0,[ones([length(IndHomo) 1]) ]);
          %%
            % running the subgrouping on individual feature
             subG_label=subG(x,y);
            % assess the fitting 
             [b0,~,~,~,stat0] = regress(y,[ones([length(IndHomo) 1]) x]);
             [b1,~,~,~,stat1] = regress(y(find(subG_label<=1.5)),...
                    [ones(size(x(find(subG_label<=1.5)))) x(find(subG_label<=1.5)) ]);
             [b2,~,~,~,stat2] = regress(y(find(subG_label>=1.5)),...
                    [ones(size(x(find(subG_label>=1.5)))) x(find(subG_label>=1.5)) ]);
   
            R(n,m,:)=[stat0(1) stat1(1) stat2(1)];
            P(n,m,:)=[stat0(3) stat1(3) stat2(3)];
            subGs(n,m,:)=subG_label;
    end
end
%%
    for n=1:Np%4005
       p_homo(n)=length(find(R(n,2:Nr+1,1)>=R(n,1,1))) / Nr;
        %R(n,1,1)
       p_heter(n)=length(find(mean(R(n,2:Nr+1,2:3),3)>=mean(R(n,1,2:3),3))) / Nr;
    end

Ind_homo = find(p_homo<=0.05);
Ind_heter = find(p_heter<=0.05);
Ind_olap = intersect(Ind_homo,Ind_heter); 

[a1,n1] = sort(R(Ind_homo,1,1),'descend');
[a2,n2] = sort(mean(R(Ind_heter,1,2:3),3),'descend');
[ao,bo] = sort(mean(R(Ind_olap,1,2:3),3),'descend');
%%
x_heter=X(IndHomo,Ind_heter);
y0=Y(IndHomo,1);
y0 = (y0-mean(y0))/std(y0);

for n=1:length(Ind_heter)
    x = x_heter(:,n);
    
    [b1,~,~,~,stat1] = regress(y0(find(subGs(Ind_heter(n),1,:)<=1.5)),...
        [ones(size(x(find(subGs(Ind_heter(n),1,:)<=1.5)))) x(find(subGs(Ind_heter(n),1,:)<=1.5)) ]);
    [b2,~,~,~,stat2] = regress(y0(find(subGs(Ind_heter(n),1,:)>=1.5)),...
        [ones(size(x(find(subGs(Ind_heter(n),1,:)>=1.5)))) x(find(subGs(Ind_heter(n),1,:)>=1.5)) ]);   
    
    for s=1:length(y0)

        k1=b1(2); 
        k2=b2(2);
    
        c1=b1(1);
        c2=b2(1);
    
       d01(n,s)=Dis_point2line(x_heter(s,n),y0(s),k1,c1);
       d02(n,s)=Dis_point2line(x_heter(s,n),y0(s),k2,c2);
    end
end
%%  whole brain level clustering 
        Dis = d01-d02;
        Mmeas=(corr(Dis(n2(1:Ntop),:),Dis(n2(1:Ntop),:)));

        DisName={'sqeuclidean'};%,correlation,''cosine,sqeuclidean
        CriName={'DaviesBouldin'};%  silhouette DaviesBouldin
        for j=1:100;
            for i=1:8
             clust(:,i,j) = kmeans(Mmeas,i,'Distance',DisName{1},'emptyaction','singleton','replicate',5);
            end
            eva = evalclusters(Mmeas,clust(:,:,j),CriName{1});
            C(j,:)=eva.CriterionValues;
            [~,indOptK_count]=min(C(j,:));
            OptK_count(j,indOptK_count)=1;
        end
[~,OptK]=min(mean(C,1)');  
%%
[~,OptClust]=min(C(:,OptK));
ClusterG=clust(:,OptK,OptClust);
[~,indSort] = sort(ClusterG,'descend');

