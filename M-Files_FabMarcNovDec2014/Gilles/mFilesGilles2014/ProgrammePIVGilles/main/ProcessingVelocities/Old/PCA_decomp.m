function [Y,var_mode]=PCA_decomp(X,var)
%======================================================================
%
% Version 1.00
% Performs a Principal Component Analysis
% 
%======================================================================
% Usage [Y,var_mode]=PCA_decomp(X,var);
% Input
%   X      :a matrix containing repeats of a measurement vector
%           each measurement vector is a column
%           X is [N x M] where M is the number of repeats
%
%   var    :variance threshold (optional)
%           the routine will output the minimum number of modes containing
%           var % of the variance in the original dat set
%           DEFAULT var=0.9.
%
% Output
%   Y        : Principal components (base of orthogonal modes in the repeats
%               of the original data X)
%   var_mode : variance contained in each mode
%
%========================================================================
%
% Update:
%       1.00    12/21/2012 Fabrice Veron
%
%========================================================================

X=X';
nvar = nargin;
if nargin==1
    var  = 0.9;
end

[m n]=size(X); %m is the number of repeats, n the length of the vector variable
h=ones(1,n);
Xmean=nanmean(X,2); %ensemble of mean-of-measurements
B=X-Xmean*h; %X with the ensemble of mean substracted out

i=isnan(B(:)); B(i)=0; %remove NaNs
C=(B*B'/n); % coviariance matrix;

[V,D] = eig(C); % V contains eigen vectors and D eignvalues on the diagonal such that X*V = V*D
lambda=diag(D); %extract eigenvalues
[lambda,ordered_indices]=sort(lambda,'descend'); %sort eigenvalues from largest to smallest
V = V(:,ordered_indices);%sort eigenvectors 

var_mode=lambda/sum(sum(lambda)); %weight for each mode
cum_norm_lambda=cumsum(var_mode); % cumulative weight 
l=min(find(cum_norm_lambda>var)); % finds lowest eigenvalues for which 90% of weight (variance) is accounted for

%makes a (reduced to the first l) sorted (ordered) subset of eigenvectors corresponding to 90% of
%variance - this is a new basis
W=V(:,1:l);

%convert the source data B to z-scores
s=sqrt(diag(C));
Z=B./(s*h);
%Z=B;
%Convert the z-score on this new (reduced) basis
W_s=conj(W');
Y=W_s*Z;
Y=Y';

F1=(W'*B);
F2=(W'*X);
F3=W_s*B;
F4=W_s*X;

var_mode=var_mode(1:l);

% Y is a reduced representation of B, it's like the "modes" of B, for example:
% Y(:,1)/var_mode(1)/100 is the same as nanmean(B);
% the modes are  Y(:,k)/(cum_norm_g(k)/100) with k=1..l
% var_mode(k)=g(k))/sum(sum(D) is the amount of variance contained in mode k
% cum_norm_g(k) is the amount of variance contained in all the modes lower than k
% here all the modes lower than k contain 90% of the variance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method 2:
% 
% %For PIV matrices where Uf contains ns snapshots of ni x nj u-velocity fields
% Uall=[reshape(Uf,ni*nj,ns)]; %Uall=[reshape(Uf,ni*nj,ns);reshape(Vf,ni*nj,ns)];
R=Uall'*Uall;         % Autocovariance matrix
[eV,D]=eig(R);        % solve: eV is eigenvectors, D is eigenvalues in diagonal matrix
[L,I]=sort(diag(D));  % sort eigenvalues in ascending order - I is sorted index vector
for i=1:length(D)
eValue(length(D)+1-i)=L(i);      % Eigenvalues sorted in descending order
eVec(:,length(D)+1-i)=eV(:,I(i)); % Eigenvectors sorted in the same order
end;
%% [L,I]=sort(diag(D),'descend'); eValue=L; eVec=eV(I);

eValue(length(eValue))=0;   % last eigenvalue should be zero
menergy=eValue/sum(eValue); % relative energy associated with mode m
l=min(find(cumsum(menergy)>var)); 
% calculate the first important modes
phi=Uall*eVec(:,1:l);  
phir=reshape(phi,ni,nj,l);

