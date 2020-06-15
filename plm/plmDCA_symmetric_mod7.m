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




function plmDCA_symmetric_mod7(fastafile,outputfile,lambda_h,lambda_J,reweighting_threshold,nr_of_cores,outputDistribution,outputMatrix)
%If should-be numericals are passed as strings, convert them.
    if (isstr(lambda_h))
        lambda_h = str2num(lambda_h);
    end
    if (isstr(lambda_J))
        lambda_J = str2num(lambda_J);
    end
    if (isstr(reweighting_threshold))
        reweighting_threshold = str2num(reweighting_threshold);
    end
    if (isstr(nr_of_cores))
        nr_of_cores = str2num(nr_of_cores);
    end

%Minimization options
    options.method='cg';	%Minimization scheme. Default: 'lbfgs', 'cg' for conjugate gradient (use 'cg' if out of RAM).
    options.progTol=1e-9;   	%Threshold for when to terminate the descent. Default: 1e-9. 

    addpath(genpath(pwd))
    
%Read inputfile (removing inserts), remove duplicate sequences, and calculate weights and B_eff.
    [N,B_with_id_seq,q,Y]=return_alignment(fastafile);
    Y=unique(Y,'rows');
    [B,N]=size(Y);
    weights = ones(B,1);
    if reweighting_threshold>0.0
        fprintf('Starting to calculate weights \n...');
        tic
        %Reweighting in MATLAB:            
        %weights = (1./(1+sum(squareform(pdist(Y,'hamm')<=reweighting_threshold))))';       
		     
        %Reweighting in C:
        Y=int32(Y);
        m=calc_inverse_weights(Y-1,reweighting_threshold);
        weights=1./m;

        fprintf('Finished calculating weights \n');
        toc
    end
    B_eff=sum(weights);
    
    [Pij_true,Pi_true] = Compute_True_Frequencies(Y,B,N,q,weights,B_eff); %computes true dist
    pseudocount_weight=0.5; %pc weight
    [Pij,Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight,N,q);
    save(outputDistribution,'Y','N','B','q','B_eff','weights','B_with_id_seq','pseudocount_weight','Pi_true','Pij_true','Pi','Pij');
    fprintf('### N = %d B_with_id_seq = %d B = %d B_eff = %.2f q = %d\n',N,B_with_id_seq,B,B_eff,q);
    
    

%Set up and run optimizer.
    field_lambda=lambda_h*B_eff;    
    coupling_lambda=lambda_J*B_eff;
    edges=[];
    for i=1:(N-1)
        for j=(i+1):N
            edges=[edges;[i,j]];
        end
    end
    Y=int32(Y);q=int32(q);edges=int32(edges);
    funObj=@(w)pseudo_likelihood_symmetric(w,Y,weights,N,q,edges,field_lambda,coupling_lambda);    
    w0=zeros(q*N+q^2*N*(N-1)/2,1);
    if nr_of_cores>1
        %matlabpool('open',nr_of_cores) 
        poolobj=parpool(nr_of_cores);
	tic
        w=minFunc(funObj,w0,options);    
        toc
        %matlabpool('close') 
	poolobj = gcp('nocreate');
	delete(poolobj);
    else
        tic
        w=minFunc(funObj,w0,options);    
        toc
    end
%Extract the estimates from w.
    htemp=reshape(w(1:q*N),N,q);    %The fields are not used in what follows.
    Jtemp=reshape(w(q*N+1:end),q,q,N*(N-1)/2);

    %convert Jtemp into N*N*q*q matrix (original gauge) 
    J=zeros(N,N,q,q);
    Counter=0;
    %top half
    for i=1:(N-1)
        for j=(i+1):N
            Counter=Counter+1;
            J(i,j,1:q,1:q)=Jtemp(1:q,1:q,Counter);
           
        end
    end
    h=htemp;
    %save(outputMatrixOriginalGauge,'J','h','N');    
    clearvars h; %important because h is htemp
    
    %bottom half
    for j=1:(N-1)
        for i=(j+1):N
            J(i,j,1:q,1:q)=reshape(J(j,i,1:q,1:q),[q,q])';
        end
    end


    
    
    %shift the fields into the Ising 
    %h_i(k)=htemp_i(k) - htemp_i(#) + sum_{j!=i} ( Jtemp_ij(k,#) - Jtemp_ij(#,#) )
    h=zeros(N,q);
    tempCouplingMatrix= zeros(q,q);
    for i=1:N
        for alpha=1:q
            h(i,alpha)=htemp(i,alpha) - mean(htemp(i,:));
            
            for j=1:N
                if j~=i
                    tempCouplingMatrix(1:q,1:q) = J(i,j,1:q,1:q);
                    h(i,alpha)=h(i,alpha)+ (mean(J(i,j,alpha,1:q)) - mean2(tempCouplingMatrix)); 
                    
                end
            end           
        end   
    end
    clearvars J; %important because J is still Jtemp (old gauge)
          
    
    

      
%A note on gauges: 
%htemp and Jtemp above satisfy the gauge
%	lambda_J*sum_s Jtemp_ij(k,s) = lambda_h*htemp_i(k)
%	lambda_J*sum_s Jtemp_ij(s,l) = lambda_h*htemp_j(l)
%	sum_s htemp_i(s) = 0.
%To obtain the full {h,J} in the Ising gauge, i.e., sum_s J_ij(k,s) = sum_s J_ij(s,l) = h_i(s) = 0, one would execute
%	J_ij(k,l)=Jtemp_ij(k,l) - Jtemp_ij(#,l) - Jtemp_ij(k,#) + Jtemp_ij(#,#)
%	h_i(k)=htemp_i(k) - htemp_i(#) + sum_{j!=i} ( Jtemp_ij(k,#) - Jtemp_ij(#,#) )
%where '#' means average. Since the fields are not used below, only the "coupling part" of this shift in implemented.

 
    %Shift the coupling estimates into the Ising gauge.
    Jising=zeros(q,q,N*(N-1)/2);
    for l=1:(N*(N-1)/2)
        Jising(:,:,l)=Jtemp(:,:,l)-repmat(mean(Jtemp(:,:,l)),q,1)-repmat(mean(Jtemp(:,:,l),2),1,q)+mean(mean(Jtemp(:,:,l)));

    end



    Counter=0;        
    J=zeros(N,N,q,q);
    for i=1:(N-1)
        for j=(i+1):N
            Counter=Counter+1;
            J(i,j,1:q,1:q)=Jising(1:q,1:q,Counter);
        end
    end

    
    %bottom half
    for j=1:(N-1)
        for i=(j+1):N
            J(i,j,1:q,1:q)=reshape(J(j,i,1:q,1:q),[q,q])';
        end
    end 
    
        
    save(outputMatrix,'J','h','N');

      
%Calculate frob. norms FN_ij.
    NORMS=zeros(N,N); 

    for i=1:(N-1)
        for j=(i+1):N
            NORMS(i,j)=norm(reshape(J(i,j,1:q,1:q),[q,q]),'fro');
            NORMS(j,i)=NORMS(i,j);

        end
    end               
       
    
%Calculate final scores, CN_ij=FN_ij-(FN_i-)(FN_-j)/(FN_--), where '-'
%denotes average.
    norm_means=mean(NORMS)*N/(N-1);
    norm_means_all=mean(mean(NORMS))*N/(N-1);
    CORRNORMS=NORMS-norm_means'*norm_means/norm_means_all;
    output=[];
    for i=1:(N-1)
        for j=(i+1):N
            output=[output;[i,j,CORRNORMS(i,j)]];
        end
    end
    dlmwrite(outputfile,output,'precision',5)
end





function A=Newmapkey(i,alpha,q) %new mapkey for mapping out q*q matrices from Nq*Nq matrix
    A = (q)*(i-1)+alpha;
end
















function [N,B,q,Y] = return_alignment(inputfile)
%Reads alignment from inputfile, removes inserts and converts into numbers.
    inputfile
    align_full = fastaread(inputfile);
    B = length(align_full);
    ind = align_full(1).Sequence ~= '.' & align_full(1).Sequence == upper( align_full(1).Sequence );
    N = sum(ind);
    Y = zeros(B,N);

    for i=1:B
        counter = 0;
        for j=1:length(ind)
            if( ind(j) )
                counter = counter + 1;
                Y(i,counter)=letter2number( align_full(i).Sequence(j) );
            end
        end
    end
    q=max(max(Y));
end

function x=letter2number(a)
    switch(a)
        % full AA alphabet
        case '-'
             x=1;
        case 'A'    
            x=2;    
        case 'C'    
            x=3;
        case 'D'
            x=4;
        case 'E'  
            x=5;
        case 'F'
            x=6;
        case 'G'  
            x=7;
        case 'H'
            x=8;
        case 'I'  
            x=9;
        case 'K'
            x=10;
        case 'L'  
            x=11;
        case 'M'
            x=12;
        case 'N'  
            x=13;
        case 'P'
            x=14;
        case 'Q'
            x=15;
        case 'R'
            x=16;
        case 'S'  
            x=17;
        case 'T'
            x=18;
        case 'V'
            x=19;
        case 'W'
            x=20;
        case 'Y'
            x=21;
        otherwise
            x=1;
    end
end





function [Pij_true,Pi_true] = Compute_True_Frequencies(align,M,N,q,W,Meff)
% Meff = B_eff above
%


% 
%     W = ones(M,1);
%     if theta>0.0
%         fprintf('Starting to calculate weights \n...');
%         tic
%         %Reweighting in MATLAB:            
%         %weights = (1./(1+sum(squareform(pdist(Y,'hamm')<=reweighting_threshold))))';       
% 		     
%         %Reweighting in C:
%         align=int32(align);
%         m=calc_inverse_weights(align-1,theta);
%         W=1./m;
% 
%         fprintf('Finished calculating weights \n');
%         toc
%     end
%     Meff=sum(W);
%     
    
    

    Pij_true = zeros(N,N,q,q);
    Pi_true = zeros(N,q);

    for j=1:M
        for i=1:N
            Pi_true(i,align(j,i)) = Pi_true(i,align(j,i)) + W(j);
        end
    end
    Pi_true = Pi_true/Meff;

    for l=1:M
        for i=1:N-1
            for j=i+1:N
                Pij_true(i,j,align(l,i),align(l,j)) = Pij_true(i,j,align(l,i),align(l,j)) + W(l);
                Pij_true(j,i,align(l,j),align(l,i)) = Pij_true(i,j,align(l,i),align(l,j));
            end
        end
    end
    Pij_true = Pij_true/Meff;

    scra = eye(q,q);
    for i=1:N
        for alpha=1:q
            for beta=1:q
                Pij_true(i,i,alpha,beta) = Pi_true(i,alpha) * scra(alpha,beta);
            end
        end
    end
end




function [Pij,Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight,N,q)
% adds pseudocount

    Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*ones(N,N,q,q);
    Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*ones(N,q);

    scra = eye(q);

    for i=1:N
        for alpha = 1:q
            for beta = 1:q
               Pij(i,i,alpha,beta) =  (1.-pseudocount_weight)*Pij_true(i,i,alpha,beta) + pseudocount_weight/q*scra(alpha,beta);
            end
        end
    end 
end






