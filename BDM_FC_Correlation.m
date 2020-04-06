load('Age_Data.mat')
Age_Data_Bin=cell(636,1);
for i=1:636
    yy=Age_Data{i};
    Age_Data_Bin{i}=Isingify2(length(yy(:,1)),498,yy);
end

t=3:258;  %time windows
FC_Corr=zeros(636,length(t));
T=ones(498,498);
T=triu(T);
T=reshape(T,[498*498,1]);
ind=find(T==0);
for i=1:length(t)
    x=randi(258-t(i)+1,1,636);  %choose random location
    for j=1:636
        %Calculate FC over a segment of fMRI data of length t(i) and
        %average over all subjects.  Repeat for the binarized (BDM) data.
        
        y=Age_Data{j}(x(j):x(j)+t(i)-1,:);
        Corr_y=corrcoef(y);
        Corr_y_reshaped=reshape(Corr_y,[498*498,1]);  
        Corr_y_reshaped=Corr_y_reshaped(ind);  %all to separate out lower triangular part
        
        y_bin=Age_Data_Bin{j}(x(j):x(j)+t(i)-2,:); %window is one shorter for BDM
        Corr_y_bin=(y_bin'*y_bin)/(t(i)-1);  %only C needs to be transposed because of shape
        Corr_y_bin_reshaped=reshape(Corr_y_bin,[498*498,1]);
        Corr_y_bin_reshaped=Corr_y_bin_reshaped(ind);
        
        % Calculate the Pearson Correlation between Corr_y and
        % Corr_y_bin
        
        R=corrcoef(Corr_y_reshaped,Corr_y_bin_reshaped);
        FC_Corr(j,i)=R(1,2);
    end
end

% Supplementary Figure ZZZ Plots FC_Corr in two different way.

% a.  The correlation between both kinds of FC was plotted vs the window
% size.  The error bars are the standard deviation across subjects (not the
% standard error of the mean).  

% b.  The correlation (with t(i) being the whole time series) between FC
% and BDM was calculated for each subject.

load('R2_DFC.mat')
scatter(AGES,FC_Corr(:,end));  %plot R2 vs Age  R values plotted not R^2.
plot(mean(FC_Corr));   %plot R2 vs TR

