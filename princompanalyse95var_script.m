% Principal Component Analysis (PCA) is a classic among the many methods of multivariate data analysis. 
% Invented in 1901 by Karl Pearson the method is mostly used today as a tool in exploratory data analysis and dimension reduction
% But also for making predictive models in machine learning.


% Reading the Input Protein Dataset-Start
% Using csvread function
% Documentation - as Per Matlab Docs
prot_raw_data=csvread('scopfinalmatlab.csv');
% Reading the Input Protein Dataset-End

% The first step of a PCA analysis is to “standardize” the data 
% The variables :
% (here the physicochemical feaures and individual amino acids/amino acids permutation 2-String) 
% are in different units and rated differently in a different scale.
% A good way to achieve this is to calculate the Z-score of our variables
% Using Z-score to normalize/standardize data.
% z score for a column of data is to substract mean/average of the column data itself
% from the individual values of features of each row and then divide the result by
% the standard deviation if the column data
% Compute zscore(X) - Start
% Documentation - as Per Matlab Docs
prot_norm_data=zscore(prot_raw_data);
% Compute Z-score(X) - End




				%Principal Component Analysis Start
						% Here we are Computing PCA using the covariance method
									% Algorithm:
										
										
										% Step1-Calculate the empirical mean
										% Step2-Calculate the deviations from the mean
										% Step3-Find the covariance matrix
										% Step4-Find the eigenvectors and eigenvalues of the covariance matrix
										% Step5-Rearrange the eigenvectors and eigenvalues Sort Descending as per index
										% Step6-Find the percentage distribution of variance per Eigenvector matrix 
										%		To find exactly how much the Eigenvector matrix Represents the whole Dataset 
										%		from where it is derived.
										% Step7-Finding the scores, i.e. get the original data 
										%		(physicochemical feaures + frequency of individual amino acids + amino acids permutation 2-String) 
										%		in term of principal components. 
										%		To do that, we just need to multiply the PCA coefficients(Eigenvector matrix) by the raw data.
										
										
										
										
										
% Aim-To find out the Principal Components


												% Step1-Calculate the empirical mean


%Compute centered variables prot_norm_data-Start
% mu is the mean/average which is a measure of central tendency. 
% uses statistics mean(X) function 
% Documentation - as Per Matlab Docs
mu = mean(prot_norm_data);
%Compute centered variables prot_norm_data-End
												
												
												
												% Step2-Calculate the deviations from the mean

% Documentation - as Per Matlab Docs
% C = bsxfun(fun,A,B) applies the element-wise binary operation specified by the function handle fun to arrays A and B.
% mu is the mean/average which is a measure of central tendency. 
% In the above code mu calculates the mean for each feature/column in the prot_norm_data matrix.
% The resultant matrix mu size is 431 X 1
% Here in the code below bsxfun function is applied to substract mean value mu 
% for each column from individual observations in the corresponding feature/column values in the prot_norm_data matrix.
% Another fact is the function bsxfun accepts arrays A and B with compatible sizes here prot_norm_data is 1059 X 431 and mu is 431 X 1
% @minus is the binary operation for substraction
prot_norm_data_dev_pca = bsxfun(@minus, prot_norm_data, mu);




												% Step3-Find the covariance matrix
% Find the Covariance
% Function used cov(X)
% Documentation - as Per Matlab Docs
covariance=cov(prot_norm_data_dev_pca);



												% Step4-Find the eigenvectors and eigenvalues of the covariance matrix
												
% Find the eigenvectors and eigenvalues of the covariance matrix covariance
% Function used eig(X)
% [V,D] = eig(A) returns matrices V and D. The columns of V present eigenvectors of A which is here covariance Matrix. 
% The diagonal matrix D contains eigenvalues
% Documentation - as Per Matlab Docs
[V,D]=eig(covariance);

												% Step5-Rearrange the eigenvectors and eigenvalues Sort Descending(eigenvalues) as per index
% Sort Descending
% Function used sort(X) as here per index i and diag(D)
% Documentation - as Per Matlab Docs	
[eigenvalues_mat,index_for_sorting]=sort(diag(D),'descend');
% Sorted Eigenvector Matrix V per index i
sorted_V=V(:,index_for_sorting);
												% Step6-Find the percentage distribution of variance per Eigenvalues matrix 
												%		To find exactly how much the Eigenvalues matrix Represents the whole Dataset 
												%		from where it is derived.
												
% Computing the percentage_distribution of the variances	
% Function used sum(x)
% Documentation - as Per Matlab Docs 											
percentage_distribution=(eigenvalues_mat./sum(eigenvalues_mat))*100;


% From PCA we want to get large variance
% Here we will aim at preserving at least 95 % variance 
% As variance only helps to determine something that varies, something that changes
% And that is where it helps to predict a class 
% Function used cumsum(X)
% Documentation - as Per Matlab Docs 
cumulative_percentage_distribution=cumsum(percentage_distribution);

% To find the first value in the 95% to 96% array of cumulative_percentage_distribution
% Handled Null in case no value exists in the afforesaid range.


% Find Logical matrix of existing values between 95 and 96 percent
least_in_95_96_range_logical=and(cumulative_percentage_distribution>95,cumulative_percentage_distribution<96);


%Find Actual values by Matrix Multiplication of cumulative_percentage_distribution and least_in_95_96_range_logical
least_in_95_96_range=cumulative_percentage_distribution(least_in_95_96_range_logical);

% Find the Size of the Matrix
select_nearest_95_variance=size(least_in_95_96_range);

% Find if the matrix contains more than one row 
if (select_nearest_95_variance(1,1)>1)
% Select the second row
    actual_variance=least_in_95_96_range(2,1);
% Build logical matrix of existing values till 1st value of 95 to 96 percent range
% As there are more than one value in between 95 and 96 percent we have to select the first value
    p_d_l_t_95_l = cumulative_percentage_distribution < actual_variance;
% Find if the matrix contains only 1 row
elseif(select_nearest_95_variance(1,1)== 1)
% Select the first row
    actual_variance=least_in_95_96_range(1,1);
% Build logical matrix of existing values till 1st value of 95 to 96 percent range 
% As there is one row means only one value between 95% and 96% so we have to select value of the row.
% So we have used less than equal to 
    p_d_l_t_95_l = cumulative_percentage_distribution <= actual_variance;
else
% This means no value exists between 95 and 96 percent so we have to go for <97 or <98 or <99
    p_d_l_t_95_l = or(cumulative_percentage_distribution < 97,cumulative_percentage_distribution < 98,cumulative_percentage_distribution < 99);
end
    
    

%Find Actual values by Matrix Multiplication of cumulative_percentage_distribution and p_d_l_t_95_l
%This will build the matrix of percentage distribution till the first value of 95% variance.
p_d_l_t_95=cumulative_percentage_distribution(p_d_l_t_95_l);

% Now this will definitely give the number of the Principal Components to select from a Data Set with 95% variance.
number_of_PCs_with_95_per_variance=size(p_d_l_t_95);
fprintf('Number of Principal Components with 95 percent variance of the selected dataset: %d \n',number_of_PCs_with_95_per_variance(1,1));

% The 95% variance is preserved in the first number_of_PCs_with_95_per_variance PCs 
% So we find the first number_of_PCs_with_95_per_variance Eigenvectors.
factors=sorted_V(:,1:number_of_PCs_with_95_per_variance);

												% Step7-Finding the scores, i.e. get the original data 
												%		(physicochemical feaures + frequency of individual amino acids + amino acids permutation 2-String) 
												%		in term of principal components. 
												%		To do that, we just need to multiply the PCA coefficients(Eigenvector matrix) by the raw data.

												
%  PCs are computed from the Eigenvectors with 95% variance.
PCA=prot_norm_data*factors;

disp('PCA written to file princompanalyse95var.csv');
dlmwrite('scopfinalmatlabpca.csv',PCA, 'precision', 9) ; % 9 significant figures.
%Principal Component Analysis End