

close all
clear all
tic %Start the timer
for mc=1:10
%-------------TYPE THERMAL SYSTEM PARAMTERS HERE--------------
%Anticipated system properties (initial guess for fitting, if you do fitting/errorbar estimation)
lambda= [243  0.1  1.32  0.028 145]; %W/m-K
C= [2.42  0.3  1.59  0.3  1.64]*1e6; %J/m^3-K
t=[58  1  135  1  1e6]*1e-9; %m 
lambda_r=[170 .1 1.32 .028 145];%isotropic layers, eta=kr/kz;
lambda_anisotropy=[0 1 0 1 0];  
lambda_ratio=lambda_r./lambda;
r=1.40e-6;

r_pump=r; %pump 1/e^2 radius, m
r_probe=r; %probe 1/e^2 radius, m
xoffset=[0 1e-6 1.5e-6];

%Standard error

lambda_error= [0 0 0 0 0]; %W/m-K
C_error= [0 0 0 0 0]*1e6; %J/m^3-K
t_error=[0 0 0 0 0]*1e-9; %m 
lambda_r_error=[0 0 0 0 0];%isotropic layers, eta=kr/kz;
r_error=0.00e-6;
r_pump_error=r_error; %pump 1/e^2 radius, m
r_probe_error=r_error; %probe 1/e^2 radius, 
phase_error=0.0;

nnodes=35;

%---2,--------- THERMAL SYSTEM PARAMTERS END HERE--------------

%choose frequencies for Sensitivity plots
freqs=logspace(log10(50e3),log10(50e6),101); %vector of frequencies (used to generate sensitivity plots)

%Choose range of frequencies to fit, Hz
freq_min=5e3;
freq_max=40e6;

%----------------------------------PROGRAM OPTIONS BEGIN--------
%Which variable(s) are you fitting for?
Xguess=[0 1 0 1 0 ;0 0 0 0 0;0 0 0 0 0 ;0 0 0 0 0;0 0 0 0 0]; %parameters to fit - this matrix has rows for lambda, C, t, eta, [r_pump r_probe dummys]. Must be the same size as input, use as many dummy 0s to pad matrix after pump/probe. Set to 1 paramenters to fit, 0 for constants.
lowbound=[ 0 0 ]; %vector of lower bound for fit parameters, same order as items above
upbound=[ 1 1];
%----------------------------------PROGRAM OPTIONS END--------


%--------------Import Data---------------------------
    files = length(xoffset);
    for K=1:files
        
        myfilename = sprintf('file%d.txt', K);
        DATAMATRIX = importdata(myfilename);
        freqs_raw=DATAMATRIX(:,1); %imported in Hertz.
        phase_raw=DATAMATRIX(:,2); %imported in degrees.
        [freqs_data(K,:),phase_data(K,:)] = extract_interior(freqs_raw',phase_raw',freq_min,freq_max);
    end
        
    
    [m,n] = size(Xguess);

        for j=1:m
            for k=1:n %iterate through the matrix
                if j==1 % this row is for lambda
                    if Xguess(j,k)==1 %if this is a variable in the fit
                        if lambda_anisotropy(k)==1
                            lambda(k)=lambda(k); %put it in place of the constant
                            lambda_r(k)=lambda(k).*lambda_ratio(k); %put it in place of the constant

                        elseif lambda_anisotropy(k)==0
                            lambda(k)=lambda(k); %put it in place of the constant
                        end
                    else
                         if lambda_anisotropy(k)==1
                            lambda(k)= lambda(k)+lambda_error(k)*randn(1,1);
                            lambda_r(k)=lambda(k).*lambda_ratio(k); %put it in place of the constant

                        elseif lambda_anisotropy(k)==0
                            lambda(k)= lambda(k)+lambda_error(k)*randn(1,1);
                        end
                    end

                elseif j==2
                         if Xguess(j,k)==1
                            C(k)=C(k);
                         else
                            C(k)= C(k)+C_error(k)*randn(1,1);
                          end
                      
                elseif j==3
                         if Xguess(j,k)==1
                            t(k)=t(k);
                         else
                            t(k)= t(k)+t_error(k)*randn(1,1);
                         end
                    
                elseif j==4
                         if Xguess(j,k)==1
                            lambda_r(k)=lambda_r(k);
                         else
                            lambda_r(k)= lambda_r(k)+lambda_r_error(k)*randn(1,1);
                         end

                elseif j==5
                         if Xguess(j,k)==1
                         if k ==1
                            r_pump(k)=r_pump(k);
                         elseif k == 2
                            r_probe(k)= r_probe(k);
                         end
                         else
                         if k ==1
                            
                            r_pump(k)= r_pump(k)+r_pump_error(k)*randn(1,1);
                            r_probe(k)= r_probe(k)+r_probe_error(k)*randn(1,1);
                     
                         end
                     end                 

                end
            end
        end
        
    phase_mixed= phase_data+phase_error*randn(size(phase_data,1),size(phase_data,2));    
        
%--------------Perform Fit (skips if no data imported)--------------------------
    fprintf('Fitting Solution:\n')
    variables = Xguess .* [lambda;C;t;lambda_r;r_pump r_probe zeros(1,size(lambda,2)-2)];
    variables = variables(find(variables));
    options=optimset('Display','Iter','TolFun',1e-6,'TolX',1e10);
    Xsol=fminsearchbnd(@(X) FDTR_FIT(X,Xguess,phase_mixed,freqs_data,lambda,C,t,lambda_r,r_pump,r_probe,nnodes,lambda_anisotropy,lambda_ratio,xoffset),variables,lowbound,upbound,options);
    fprintf('Data fit completed\n')
    Xsol';
    for K=1:length(xoffset)
        [~, phase_model(K,:), deltaT_model(K,:)] = FDTR_FIT(Xsol,Xguess,phase_mixed,freqs_data,lambda,C,t,lambda_r,r_pump,r_probe,nnodes,lambda_anisotropy,lambda_ratio,xoffset(K));
   
    end

     for j=1:size(variables,1)
        Xsol_up=Xsol';
        Xsol_up(j)=Xsol_up(j)*1.01; %calculate model with 1% variation of each variable
        x=[];
        for k=1:length(xoffset)
            [~, temp, ~] = FDTR_FIT(Xsol_up,Xguess,phase_data,freqs_data,lambda,C,t,lambda_r,r_pump,r_probe,nnodes,lambda_anisotropy,lambda_ratio,xoffset(k));
            
                x=cat(1, x, temp');
      
        end
        up(:,j)=x;
    end
    
    dof=size(freqs_data,2)*length(xoffset)-size(variables,1); % degrees of freedom
    s=0;
    phase_model_J=[];
    
    for k=1:length(xoffset)
        s=s+sum((phase_data(k,:)-phase_model(k,:)).^2);
        sdr=sqrt(s/dof); %standard deviation of the residulas
        phase_model_J=cat(1,phase_model_J,phase_model(k,:)');
		
    end
    
    % Jacobian matrix
    for j=1:size(variables,1)
        J(:,j)=(up(:,j)-phase_model_J)./(Xsol_up(j)*0.01);
    end
    
    % standard error
    CX=pinv(J'*J);
    sigma=sdr^2*CX;
    se=sqrt(diag(sigma))';
    Xsolution(mc,:)=Xsol;
    serror(mc,:)=se;
end  

    figure(1)
    semilogx(freqs_data',phase_model','*','MarkerSize',8)
    hold on
    semilogx(freqs_data',phase_mixed','o','MarkerSize',8)
    set(gca,'FontSize',16)
    xlabel('Frequency (Hz)','Fontsize',16)
    ylabel('Phase (deg)','FontSize',16)

toc

fprintf('Program Completed\n')

[x,y]=size(Xsolution);
for h=1:y
    figure(h+1)
    hist(Xsolution(:,h),20)
end
serror=serror';
Xsolution=Xsolution';
for i=1:length(Xsolution)
    w(:,i)=1./(serror(:,i).*serror(:,i));  
end
n=w.*Xsolution;
nu=sum(n');
den=sum(w');
mean=nu./den
Uncertaininty_on_mean=1./sqrt(den);

