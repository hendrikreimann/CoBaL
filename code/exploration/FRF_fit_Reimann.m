%	FRF_fit_Reimann.m		27-May-2021%function [Fit,fval]=FRF_fit_Reimann(F,FRF,J,mgh,Nfit)%% Inputs: F: vector of FRF frequencies (Hz)%         FRF: vector of FRF values (complex values)%         J: Subject moment of inertia about ankle joint%         mgh: Subject mass (excluding feet) times gravity constant times%                CoM height above ankle joint%         Nfit: number of optimization fits to perform using different%                initial values of parameters%% Output: Fit: Structure variable holding best fit parameter values%         fval: Fit error value of best fit%ferr=@(P_fit_i)FRF_fit_err_Reimann(P_fit_i,F,FRF,J,mgh); % function to calculate fit error% random_lb and random_ub are used to calculate initial parameter values% parameters are [gn(sensory weight) kp(proportional gain) kd(derivative gain) kf(force feedback) td(time delay)]random_lb=[0.1 1.2*mgh 0.3*mgh 0.00005 0.1]; % lower boundrandom_ub=[0.9 1.6*mgh 0.6*mgh 0.0002  0.2]; % upper bound% random_lb=[0.6 1.5*mgh   0.5*mgh 0.0001 0.12];% random_ub=[0.6 1.5*mgh   0.5*mgh 0.0001 0.12];%% fit_lb and fit_ub are the constrained absolute search limits used by the% optimization functionfit_lb=[0.001 mgh  0.1*mgh 0.00000001  0.02];	% fit lower boundsfit_ub=[101.1  4*mgh 2*mgh 4*pi/180   0.35];	% fit upper boundsA=[];B=[];Aeq=[];Beq=[];options=optimset('fmincon');options=optimset(options,'MaxFunEvals',10000,'TolFun',1e-10,'MaxIter',1000,'Display','off');% options=optimset(options,'MaxFunEvals',10000,'TolFun',1e-10,'MaxIter',1000);%options=optimset('MaxFunEvals',10000,'TolFun',1e-10,'MaxIter',1000,'Display','off');%% Perform Nfit fits and save the one with the lowest error%fval=100;   % start with extremely large error valuefor ii=1:Nfit    P_fit_i=random_lb+rand(size(random_lb)).*(random_ub-random_lb); % initial values    [P_fit_j,fval_j]=fmincon(ferr,P_fit_i,A,B,Aeq,Beq,fit_lb,fit_ub,[],options); % Optimization Toolbox function fmincon    if fval_j<fval        P_fit=P_fit_j;        fval=fval_j;    endendFit.gn=P_fit(1);Fit.kp=P_fit(2); % in radian unitsFit.kd=P_fit(3);Fit.kf=P_fit(4);Fit.td=P_fit(5);