function [Z, phase_model, deltaT_model]=FDTR_FIT(X,XGuess,phase_mixed,freqs,lambda,C,t,lambda_r,r_pump,r_probe,nnodes,lambda_anisotropy,lambda_ratio,xoffset)

%Define the variables to be fit
%Look at XGuess matrix, and for elements = 1 substitute the current fit
%vector into the list of arguments to be modeled. i.e. if XGuess shows that
%we are fitting for lambda(1), then substitute the current value in X th
%represents lambda(1) into the lambda vecotr before calling FDTR_REFL
[m,n] = size(XGuess);
count = 1; %counts parameter number in vector X being evaluated
for j=1:m
    for k=1:n %iterate through the matrix
        if j==1 % this row is for lambda
            if XGuess(j,k)==1 %if this is a variable in the fit
				if lambda_anisotropy(k)==1
					lambda(k)=X(count); %put it in place of the constant
					lambda_r(k)=lambda(k).*lambda_ratio(k); %put it in place of the constant
					count = count + 1; %keep track of the parameter number
				elseif lambda_anisotropy(k)==0
					lambda(k)=X(count); %put it in place of the constant
					count = count + 1; %keep track of the parameter number
				end
            end
        elseif j==2
            if XGuess(j,k)==1
                C(k)=X(count);
                count = count + 1;
            end
        elseif j==3
            if XGuess(j,k)==1
                t(k)=X(count);
                count = count + 1;
            end
        elseif j==4
            if XGuess(j,k)==1
                lambda_r(k)=X(count);
                count = count + 1;
            end
        elseif j==5
            if XGuess(j,k)==1
                if k == 1
                    r_pump=X(count);
                elseif k == 2
                    r_probe=X(count);
                end
                count = count + 1;
       
            end
        end
    end
end

% lambda(3)=X(1);
% lambda(2)=X(2);




Z=0;
for K =1:length(xoffset)
    [deltaT_model(K,:),phase_model(K,:)]=FDTR_REFL(freqs(K,:),lambda,C,t,lambda_r,r_pump,r_probe,nnodes,lambda_anisotropy,lambda_ratio,xoffset(K));
    res=(phase_model(K,:)-phase_mixed(K,:)).^2;
    Z=Z+sum(res);
end