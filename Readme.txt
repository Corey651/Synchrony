This code accompanies the publication "Metabolism modulates network synchrony in the aging brain".

1.  Diet_Data has the pre-formated matlab data (Diet experiment) in the form of 
	a 12x2 cell array, corresponding to our 12 subjects and our two diets
	(glycolytic and ketogenic).  Each entry of the cell is a data matrix 
	containing the fMRI data for each subject:
	 740 (TR) x 498 (Regions)

2.  Age_Data has the pre-formated matlab data (Camcan) in a 636x1
	cell array, corresponding to the 636 subjects from the publicly available
	Camcan dataset.  Each entry of the cell matrix is a data matrix like above:
	258 (TR) x 498 (Regions).  The ages for each of these subjects is contained in
	Sub_Ages.mat

*** To produce Age_Data, run Combine_Age_Data.m to combine the four parts of age data
	together into a single file.

3.  Sub_Ages has the ages of each subject (in Camcan, 636 entries)

4.  Flip_Ind (or more specifically the first 16 entries) contains the 
        regions whose signs we flip.

5.  Find_Clusters.m computes Flip_Ind

6.  Fit_Lambda_Script.m computes all lambda values and associated statistics.

7.  lamage, lamket, lamglu contain (in order) the best fit (raw) lambda 
	for the lifespan data and diet (ketogenic and glycolytic).  The rescaled Lambda
	often used in the text is Lambda=(lambda-lam_crit)/lam_crit (see text).