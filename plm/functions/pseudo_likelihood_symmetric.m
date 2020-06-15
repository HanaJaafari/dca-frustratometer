% Copyright 2014 - by Magnus Ekeberg (magnus.ekeberg@gmail.com)
% All rights reserved
% 
% Permission is granted for anyone to copy, use, or modify this
% software for any uncommercial purposes, provided this copyright 
% notice is retained, and note is made of any changes that have 
% been made. This software is distributed without any warranty, 
% express or implied. In no event shall the author or contributors be 
% liable for any damage arising out of the use of this software.
% 
% The publication of research using this software, modified or not, must include 
% appropriate citations to:
%
% 	M. Ekeberg, C. Lövkvist, Y. Lan, M. Weigt, E. Aurell, Improved contact
% 	prediction in proteins: Using pseudolikelihoods to infer Potts models, Phys. Rev. E 87, 012707 (2013) 
%
%	M. Ekeberg, T. Hartonen, E. Aurell, Fast pseudolikelihood
%	maximization for direct-coupling analysis of protein structure
%	from many homologous amino-acid sequences, arXiv:1401.4832
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [pseudoNLL,g] = pseudo_likelihood_symmetric(w,Y,weights,N,q,edges,lambdah,lambdaJ)
%This function evaluates the full pseudolikelihood as a sum of marginal likelihoods g_r.
%Unlike in the paper's formulation, the regularization terms are here included in g_r.

%Prepare inputs to each g_r; g_r does not take all parameters h and J, but only h_r and J_r (={J_ri}_[i!=r]).
h=reshape(w(1:q*N),N,q); J=reshape(w(q*N+1:end),q,q,N*(N-1)/2);
J_rs=zeros(q^2*(N-1),N); %Matrix in which to store couplings for input to g_r (with column r containing J_r).
for r=1:N
   ls1=(edges(:,2)==r);ls2=edges(:,1)==r; %Extract all pairs that have the position index r in them.
   J_rs(1:q^2*(r-1),r)=reshape(permute(J(:,:,ls1),[2 1 3]),q^2*(r-1),1); 
   J_rs(q^2*(r-1)+1:end,r)=reshape(J(:,:,ls2),q^2*(N-r),1);
%The reason for the permutations above: 
%For J_ij as stored in the variable J, i is always <j. But, in g_r the first AA index is always associated with position r, as in J_ri(sigma_r,sigma_i). So, if r>i sigma_r and sigma_i must be swapped.
end

%Evaluate all g_r.
f=zeros(1,N);
grad=zeros(q+q^2*(N-1),N);
parfor r=1:N
    [f(r) grad(:,r)]=g_r([h(r,:)';J_rs(:,r)],Y,weights,N,q,lambdah,lambdaJ/2,r);	%Since a regularization term for J_ij is included in both g_i and g_j, we divide the coupling-regularization strength by 2 to compensate for double counting.
end

%Assemble the full gradient.
grad1=grad(1:q,:)';
grad_reshaped=reshape(grad(q+1:end,:),q,q,N-1,N);
grad2=zeros(q,q,N*(N-1)/2);
l=1;
for i=1:(N-1)
    for j=(i+1):N
        grad2(:,:,l)=grad_reshaped(:,:,j-1,i)+grad_reshaped(:,:,i,j)'; %Contributions to the derivative w.r.t. J_ij come from g_i and g_j.
        l=l+1;
    end
end

pseudoNLL=sum(f); g=[grad1(:);grad2(:)];
