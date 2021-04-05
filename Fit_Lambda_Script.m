%Here we start with all of the data and fit all lambda values
%Pilot Diet data, N=12 subjects, length 740 TR
%Camcan Lifespan data, N=636, length 258 TR

%Here we do calculations with the raw lambda value (as opposed to
%Lambda=(lambda-lambda_c)/lambda_c).


%% Binarizing the Data
load('Diet_Data.mat');                       %load variable Diet_Data

load('Flip_Ind.mat');                        %load variable Flip_Ind

%The first 16 entries here are the ROI we decided to flip the signs of for 

Flip=ones(1,498);                                 %1 if we don't flip, -1 if we flip.  See Methods and XYZ.m
for i=1:15
    Flip(Flip_Ind(i))=-1;
end


Bglu=zeros(739,12);                              
Bket=zeros(739,12);
for i=1:12
    yy=Diet_Data{i,1};
    yyn=Isingify2(length(yy(:,1)),498,yy);         %Isingify2.m binarizes input data using BDM
    Bglu(:,i)=sum(Flip.*yyn,2)/498;
    
    yy=Diet_Data{i,2};
    yyn=Isingify2(length(yy(:,1)),498,yy);
    Bket(:,i)=sum(Flip.*yyn,2)/498;    
end
load('Age_Data.mat');                              %Loads Age_Data
BAge=zeros(257,636);
for i=1:636
    yy=Age_Data{i};
    yyn=Isingify2(length(yy(:,1)),498,yy);
    BAge(:,i)=sum(Flip.*yyn,2)/498;
end
%%

%Here we fit Lambda (Lam) using the probability distribution Eq 5

lamglu=zeros(1,12);            
lamket=zeros(1,12);
lamage=zeros(1,636);


     %This finds the unique value of lambda that best fits the sufficient
     %statistic, mean(s^2), for each subject
     
for t=1:12
    m2=mean(Bglu(:,t).^2);                                      
    k=0:1:498;
    v=(2*k-498)/498;  %all values of s
    vv=v.^2; 
    nck=zeros(1,499);
    for i=1:499
        nck(i)=nchoosek(498,k(i));
    end
    f1=@(lambda) sum((vv).*nck.*exp(lambda.*vv*498^2));
    f2=@(lambda) sum(nck.*exp(lambda.*vv*498^2));               %ML estimation
    f3=@(lambda) f1(lambda)./f2(lambda)-m2;
    options = optimset('TolX',1*10^-10);
    lamglu(t)=fzero(f3,.00001,options);          
end


for t=1:12
    m2=mean(Bket(:,t).^2);
    k=0:1:498;
    v=(2*k-498)/498;
    vv=v.^2;
    nck=zeros(1,499);
    for i=1:499
        nck(i)=nchoosek(498,k(i));
    end
    f1=@(lambda) sum((vv).*nck.*exp(lambda.*vv*498^2));
    f2=@(lambda) sum(nck.*exp(lambda.*vv*498^2));
    f3=@(lambda) f1(lambda)./f2(lambda)-m2;
    options = optimset('TolX',1*10^-10);
    lamket(t)=fzero(f3,.00001,options);
end


for t=1:636
    m2=mean(BAge(:,t).^2);
    k=0:1:498;
    v=(2*k-498)/498;
    vv=v.^2;
    nck=zeros(1,499);
    for i=1:499
         nck(i)=nchoosek(498,k(i));
    end
    f1=@(lambda) sum((vv).*nck.*exp(lambda.*vv*498^2));
    f2=@(lambda) sum(nck.*exp(lambda.*vv*498^2));
    f3=@(lambda) f1(lambda)./f2(lambda)-m2;
    options = optimset('TolX',1*10^-10);
    lamage(t)=fzero(f3,.00001,options);
end

%Outputs up to this point are ML lambda values for all subjects


    %Statistics based on ML estimates
 
load('Sub_Ages.mat');                              %Loads Sub_Ages

pdiet=signrank(lamglu,lamket,'tail','left');       %Wilcoxon Sign-rank for diets

[~,page]=corr(Sub_Ages,lamage','Type','Spearman'); %Spearman-Rank Test for age

%statistics were calculated using the raw lamda.  Using
%Lambda=(lambda-lambda_c)/(lambda_c) instead does not change the pvalues or
%test outcomes.

%% Find Errorbars for Each Lambda fit


x=.0006:.00001:.0016;                           %Essentially our uniform prior distribution,  other lambda values are unphysical

%Calculating the log-likelihood ratio of the ML fit from the data for each
%subject to all of the lambda values in x

LLglur=zeros(12,length(x));
LLketr=zeros(12,length(x));
LLager=zeros(636,length(x));

for i=1:636
    va=BAge(:,i);
    Qa=sum(va.^2);
    for j=1:length(x)
        LLager(i,j)=(x(j)-lamage(i))*Qa*498^2-length(va)*(log(f2(x(j)))-log(f2(lamage(i))));
    end
end


for i=1:12
    vglu=Bglu(:,i);
    vket=Bket(:,i);
    Qsglu=sum(vglu.^2);
    Qsket=sum(vket.^2);
    for j=1:length(x)
        LLglur(i,j)=(x(j)-lamglu(i))*Qsglu*498^2-length(vglu)*(log(f2(x(j)))-log(f2(lamglu(i))));
        LLketr(i,j)=(x(j)-lamket(i))*Qsket*498^2-length(vket)*(log(f2(x(j)))-log(f2(lamket(i))));
    end
end

%Calculating the actual errorbars

Errtopg=zeros(1,12);
Errbotg=zeros(1,12);

Errtopk=zeros(1,12);
Errbotk=zeros(1,12);

Errtopage=zeros(1,636);
Errbotage=zeros(1,636);

thr=log(.01);             %LL ratio is .01, that means that the ends of our confidence interval are predicted 
                          %to be 100-fold less likely than the ML estimate.
                
%for each subject, finds the index (for the top and the bottom) of x, and corresponding values of lambda
%where the likelihood ratio is .01.  This is done by looping through x and
%stopping at the lowest and highest values in x that satify LL(x(i))>.01

for i=1:636
    xx=x-lamage(i);
    qp=find(xx>0);
    qm=find(xx<0);
    ind=0;
    ird=1;
    while ind==0
        if LLager(i,qp(ird))<=thr  
            Errtopage(i)=x(qp(ird));
            ind=1;
        elseif ird>length(qp)
            ind=1;
        else
            ird=ird+1;
        end
    end
    ind=0;
    ird=1;
    while ind==0
        if LLager(i,qm(length(qm)-ird+1))<=thr
            Errbotage(i)=x(qm(length(qm)-ird+1));
            ind=1;
        elseif ird>length(qm)
            ind=1;
        else
            ird=ird+1;
        end
    end
end
            
for i=1:12
    xx=x-lamglu(i);
    qp=find(xx>0);
    qm=find(xx<0);
    ind=0;
    ird=1;
    while ind==0
        if LLglur(i,qp(ird))<=thr
            Errtopg(i)=x(qp(ird));
            ind=1;
        elseif ird>length(qp)
            ind=1;
        else
            ird=ird+1;
        end
    end
    ind=0;
    ird=1;
    while ind==0
        if LLglur(i,qm(length(qm)-ird+1))<=thr
            Errbotg(i)=x(qm(length(qm)-ird+1));
            ind=1;
        elseif ird>length(qm)
            ind=1;
        else
            ird=ird+1;
        end
    end
end

for i=1:12
    xx=x-lamket(i);
    qp=find(xx>0);
    qm=find(xx<0);
    ind=0;
    ird=1;
    while ind==0
        if LLketr(i,qp(ird))<=thr
            Errtopk(i)=x(qp(ird));
            ind=1;
        elseif ird>length(qp)
            ind=1;
        else
            ird=ird+1;
        end
    end
    ind=0;
    ird=1;
    while ind==0
        if LLketr(i,qm(length(qm)-ird+1))<=thr
            Errbotk(i)=x(qm(length(qm)-ird+1));
            ind=1;
        elseif ird>length(qm)
            ind=1;
        else
            ird=ird+1;
        end
    end
end


%output Error bars for all Lambda fits


%% Find probability distribution of Ising model 

k=0:1:498;
v=(2*k-498)/498;
vv=v.^2;
nck=zeros(1,499);
for i=1:499
    nck(i)=nchoosek(498,k(i));
end
x=.0007:.00001:.0012;
f1=@(lambda) sum(nck.*exp(lambda.*vv*498^2));
y=zeros(499,length(x));
for j=1:499
    for i=1:length(x)
        y(j,i)=x(i)*((2*(j-1)-498)^2)-log(f1(x(i)));
        y(j,i)=y(j,i)+log(nchoosek(498,498+1-j));
    end
end
P=sum(exp(y),1);
PP=log(P);
Prob_Dist=y-PP;  %normalization correction just in case
Prob_Dist=exp(Prob_Dist);   %This gives the probabilities of each synchrony for different values of lambda

% Calculating the 2nd and 4th moments
 
    % from data- these were all concatenated in Supplementary Fig. ABC 
    
K2g=mean(Bglu.^2,1);
K4g=var(Bglu.^2,1,1);

K2k=mean(Bket.^2,1);
K4k=var(Bket.^2,1,1);

K2a=mean(BAge.^2,1);
K4a=var(BAge.^2,1,1);

    %from the model

K2mod=sum(vv'.*Prob_Dist,1);
K4mod=sum(vv'.^2.*Prob_Dist,1)-K2mod.^2;

%% Calculating the means and error bars from binning each subject.

   % For the plots we used the rescaled Lambda, so we will compute that
   % here
   
lam_crit=1.004*10^(-3);
   
Lamglu=(lamglu-lam_crit)/lam_crit;          %Capital L refers to the rescaled values
Lamket=(lamket-lam_crit)/lam_crit;
Lamage=(lamage-lam_crit)/lam_crit;


   % Age groups
   

Edges=[18,32,46,60,74,88];  %The plotted age groups (5 total) are the midpoints between these edges 

Q25e=zeros(1,5);
Q75e=zeros(1,5);
Qmed=zeros(1,5);

for i=1:5
    q=find(Sub_Ages>Edges(i)&Sub_Ages<Edges(i+1));          %calculates the medians and quartiles (error bars) of each age group 
                                                            % Supplementary Fig Z
    sss=Lamage(q);
    Qmed(i)=median(sss);
    Q25e(i)=Qmed(i)-quantile(sss,.25);
    Q75e(i)=quantile(sss,.75)-Qmed(i);
end

    % Diet
    
    %Here we worked with difLam=Lamket-Lamglu
    
difLam=Lamket-Lamglu;

    %Next we have to propagate the error bars (from our posterior
    %distribution to find the error bars of difLam for each subject
    
    % For simplicity we just approximate the error as being symmetric.
  
Errg=(Errtopg-Errbotg)/2;
Errk=(Errtopk-Errbotk)/2;
difErr=sqrt(Errg.^2+Errk.^2)/(lam_crit);            %standard error propagation for differences
                                                    %rescaling was because of rescaling lambda->Lambda
                                                    
%difErr gives the errorbars in Fig 5, where the
%subjects were ordered by difLam.

%% Find Int and Seg and find their ratios vs Lambda


All_Age_data=zeros(257*636,498);
for i=1:636
    yy=Age_Data{i};
    yyn=Isingify2(length(yy(:,1)),498,yy);
    All_Age_data(257*(i-1)+1:257*(i-1)+257,:)=Flip.*yyn;
end
BAge_reshaped=reshape(BAge,[257*636,1]);

%Find Int and Seg and normalize them
Ind_Seg=All_Age_data((abs(BAge_reshaped)==0),:);
Seg=corrcoef(Ind_Seg);
Seg=Seg-eye(498,498);
Seg=Seg/sqrt(sum(sum(Seg.^2)));  
Ind_Int=All_Age_data((abs(BAge_reshaped)>=1/2),:);
Int=corrcoef(Ind_Int);
Int=Int-eye(498,498);
Int=Int/sqrt(sum(sum(Int.^2)));

%Find P_Seg and P_Int for each subject

pseg=zeros(1,636);
pint=zeros(1,636);

cov_Seg_Int=sum(sum(Seg.*Int));
for i=1:636
    yy=Age_Data{i};
    yyn=Flip.*Isingify2(length(yy(:,1)),498,yy);
    Ma=cov(yyn);
    Ma=Ma-diag(diag(Ma));
    Ma=Ma/sqrt((sum(sum(Ma.^2))));
    pinti=sum(sum(Ma.*Int));
    psegi=sum(sum(Ma.*Seg));
    pint(i)=(pinti-psegi*cov_Seg_Int)/(1-cov_Seg_Int^2);
    pseg(i)=(psegi-pinti*cov_Seg_Int)/(1-cov_Seg_Int^2);
end

pintd=zeros(1,24);
psegd=zeros(1,24);
for i=1:12
    yy=Diet_Data{i,1};
    yyn=Flip.*Isingify2(length(yy(:,1)),498,yy);
    Ma=cov(yyn);
    Ma=Ma-diag(diag(Ma));
    Ma=Ma/sqrt(sum(sum(Ma.^2)));
    pinti=sum(sum(Ma.*Int));
    psegi=sum(sum(Ma.*Seg));
    pintd(i)=(pinti-psegi*cov_Seg_Int)/(1-cov_Seg_Int^2);
    psegd(i)=(psegi-pinti*cov_Seg_Int)/(1-cov_Seg_Int^2);
    yy=Diet_Data{i,2};
    yyn=Flip.*Isingify2(length(yy(:,1)),498,yy);
    Ma=cov(yyn);
    Ma=Ma-diag(diag(Ma));
    Ma=Ma/sqrt(sum(sum(Ma.^2)));
    pinti=sum(sum(Ma.*Int));
    psegi=sum(sum(Ma.*Seg));
    pintd(i+12)=(pinti-psegi*cov_Seg_Int)/(1-cov_Seg_Int^2);
    psegd(i+12)=(psegi-pinti*cov_Seg_Int)/(1-cov_Seg_Int^2);
end

Lam_All=[Lamage,Lamglu,Lamket];
pint=[pint,pintd];
pseg=[pseg,psegd];

pint=pint./(pint+pseg);
pseg=1-pint;

%% Plot typical fits (for lambda above, below, and at the critical point)
load('lamket.mat')
load('lamglu.mat')
load('lamage.mat')


% Age 32

figure

annotation('textarrow','Position',[.0675,0.085,0,.15],'String','   P(s)','HeadWidth',5,'HeadLength',5,'HorizontalAlignment','left','VerticalAlignment','cap','TextRotation',90,'Textmargin',10)
annotation('doublearrow','Position',[.0250,0.085,.075,0],'Head1Width',5,'Head2Width',5,'Head1Length',5,'Head2Length',5)
annotation('textbox',[.0475,0.055,0,.05],'String','s','FitBoxToText','on','EdgeColor','none')

Lamplot=lamage(120);
Q=120;
pfit=zeros(1,499);
k=0:1:498;
v=(2*k-498)/498;  %all values of s
vv=v.^2; 
nck=zeros(1,499);
for i=1:499
    nck(i)=nchoosek(498,k(i));
    pfit(i)=nck(i)*exp(Lamplot*vv(i)*498^2);
end
countsp=round(pfit);
[countsd,Edges]=histcounts(BAge(:,Q));
countspnew=zeros(1,length(Edges)-1);
for i=1:length(countspnew)
    k=ceil(Edges(i)*498);
    if mod(k,2)==1
        k=k+1;
    end
    while k<= floor(Edges(i+1)*498)
        countspnew(i)=countspnew(i)+countsp((k+500)/2);
        k=k+2;
    end
end

Edges1=Edges;
countspnew1=countspnew;

subplot(2,3,1)
hold on
histogram('BinEdges',Edges1,'BinCounts',countspnew1,'normalization','probability','DisplayStyle','stairs');
histogram(BAge(:,Q),'normalization','probability','DisplayStyle','stairs')
xlabel('Age 32')

% Age 48

Lamplot=lamage(250);
Q=250;
pfit=zeros(1,499);
k=0:1:498;
v=(2*k-498)/498;  %all values of s
vv=v.^2; 
nck=zeros(1,499);
for i=1:499
    nck(i)=nchoosek(498,k(i));
    pfit(i)=nck(i)*exp(Lamplot*vv(i)*498^2);
end
countsp=round(pfit);
[countsd,Edges]=histcounts(BAge(:,Q));
countspnew=zeros(1,length(Edges)-1);
for i=1:length(countspnew)
    k=ceil(Edges(i)*498);
    if mod(k,2)==1
        k=k+1;
    end
    while k<= floor(Edges(i+1)*498)
        countspnew(i)=countspnew(i)+countsp((k+500)/2);
        k=k+2;
    end
end

Edges2=Edges;
countspnew2=countspnew;

subplot(2,3,2)
hold on
histogram('BinEdges',Edges2,'BinCounts',countspnew2,'normalization','probability','DisplayStyle','stairs');
histogram(BAge(:,Q),'normalization','probability','DisplayStyle','stairs')
xlabel('Age 48')


% Age 63
Lamplot=lamage(450);
Q=450;
pfit=zeros(1,499);
k=0:1:498;
v=(2*k-498)/498;  %all values of s
vv=v.^2; 
nck=zeros(1,499);
for i=1:499
    nck(i)=nchoosek(498,k(i));
    pfit(i)=nck(i)*exp(Lamplot*vv(i)*498^2);
end
countsp=round(pfit);
[countsd,Edges]=histcounts(BAge(:,Q));
countspnew=zeros(1,length(Edges)-1);
for i=1:length(countspnew)
    k=ceil(Edges(i)*498);
    if mod(k,2)==1
        k=k+1;
    end
    while k<= floor(Edges(i+1)*498)
        countspnew(i)=countspnew(i)+countsp((k+500)/2);
        k=k+2;
    end
end

Edges3=Edges;
countspnew3=countspnew;

subplot(2,3,4)
hold on
histogram('BinEdges',Edges3,'BinCounts',countspnew3,'normalization','probability','DisplayStyle','stairs');
histogram(BAge(:,Q),'normalization','probability','DisplayStyle','stairs')
xlabel('Age 63')



% Age 81

Lamplot=lamage(600); 
Q=600;
pfit=zeros(1,499);
k=0:1:498;
v=(2*k-498)/498;  %all values of s
vv=v.^2; 
nck=zeros(1,499);
for i=1:499
    nck(i)=nchoosek(498,k(i));
    pfit(i)=nck(i)*exp(Lamplot*vv(i)*498^2);
end
countsp=round(pfit);
[countsd,Edges]=histcounts(BAge(:,Q));
countspnew=zeros(1,length(Edges)-1);
for i=1:length(countspnew)
    k=ceil(Edges(i)*498);
    if mod(k,2)==1
        k=k+1;
    end
    while k<= floor(Edges(i+1)*498)
        countspnew(i)=countspnew(i)+countsp((k+500)/2);
        k=k+2;
    end
end

Edges4=Edges;
countspnew4=countspnew;

subplot(2,3,5)
hold on
histogram('BinEdges',Edges4,'BinCounts',countspnew4,'normalization','probability','DisplayStyle','stairs');
histogram(BAge(:,Q),'normalization','probability','DisplayStyle','stairs')
xlabel('Age 81')


% Glucose

Q=10;

Lamplot=lamglu(Q);
pfit=zeros(1,499);

k=0:1:498;
v=(2*k-498)/498;  %all values of s
vv=v.^2; 
nck=zeros(1,499);
for i=1:499
    nck(i)=nchoosek(498,k(i));
    pfit(i)=nck(i)*exp(Lamplot*vv(i)*498^2);
end
countsp=round(pfit);
[countsd,Edges]=histcounts(Bglu(:,Q));
countspnew=zeros(1,length(Edges)-1);
for i=1:length(countspnew)
    k=ceil(Edges(i)*498);
    if mod(k,2)==1
        k=k+1;
    end
    while k<= floor(Edges(i+1)*498)
        countspnew(i)=countspnew(i)+countsp((k+500)/2);
        k=k+2;
    end
end

Edges5=Edges;
countspnew5=countspnew;

subplot(2,3,3)
hold on
histogram('BinEdges',Edges5,'BinCounts',countspnew5,'normalization','probability','DisplayStyle','stairs');
histogram(Bglu(:,Q),'normalization','probability','DisplayStyle','stairs')
xlabel('Glucose')

% Ketones

Q=10;

Lamplot=lamket(Q);
pfit=zeros(1,499);

k=0:1:498;
v=(2*k-498)/498;  %all values of s
vv=v.^2; 
nck=zeros(1,499);
for i=1:499
    nck(i)=nchoosek(498,k(i));
    pfit(i)=nck(i)*exp(Lamplot*vv(i)*498^2);
end
countsp=round(pfit);
[countsd,Edges]=histcounts(Bket(:,Q));
countspnew=zeros(1,length(Edges)-1);
for i=1:length(countspnew)
    k=ceil(Edges(i)*498);
    if mod(k,2)==1
        k=k+1;
    end
    while k<= floor(Edges(i+1)*498)
        countspnew(i)=countspnew(i)+countsp((k+500)/2);
        k=k+2;
    end
end

Edges6=Edges;
countspnew6=countspnew;

subplot(2,3,6)
hold on
histogram('BinEdges',Edges6,'BinCounts',countspnew6,'normalization','probability','DisplayStyle','stairs');
histogram(Bket(:,Q),'normalization','probability','DisplayStyle','stairs')
xlabel('Ketones')

