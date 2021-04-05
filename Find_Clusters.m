%% Finds the Postive and Negative clusters

% Here we determine which regions to flip the sign of.  We do this by:

% First: Binarizing our data

load('Age_Data.mat');    %variable Age_Data
Bin_Age=cell(636,1);
BAge=zeros(257,636);
for i=1:636
    yy=Age_Data{i};
    Bin_Age{i}=Isingify2(length(yy(:,1)),498,yy);
end

% Second: Calculating the age-averaged correlations for all regions

CorrMat=zeros(498,498);
for i=1:636
    CorrMat=CorrMat+corrcoef(Bin_Age{i});
end

CorrMat=CorrMat/636;
% Third: Calculating the average correlation per region

TotCor=sum(CorrMat);
[TotCor_Ordered,Flip_Ind]=sort(TotCor);

        % Plotting TotCor shows that there are several obvious outliers
        % with very negative total correlations.  These are the regions
        % that we flipped. (See Supplementary Fig YYY).
        
        % We flipped the first 16 regions; representing three distinct
        % functional domains; a mixture subcortical regions (thalamus,
        % hippocampus, caudate nucleus), the supplementary motor area, and
        % the lingual gyrus.
        
% Saved in Flip_Ind.mat
        