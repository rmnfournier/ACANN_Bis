%% Maximum Entropy Implementation 
clear all;close all;
%% Parameters
Beta=1;
maxiter=10000;
singular_value_tol=1e-7;
mean_data=[];
var_data=[];
data = [];
tau=[];
folder=".";
for i=0:31
    filename=folder+"/simulation"+int2str(i)+".csv";
    tmp=importdata(filename);
    mean_data=[mean_data, tmp(:,2)];
    if(length(tmp(1,:))<3)
           var_data=[var_data,0];
    else
        var_data=[var_data,tmp(:,3)];
    end
    tau=[tau,(tmp(:,1))];
    filename=folder+"/correlations"+int2str(i)+".csv";
    data=[data;importdata(filename)];
end
var_data=var_data/mean_data(1);
mean_data=mean_data/mean_data(1);

nb_alphas=200;
alpha=1;
step_alpha=1.3;
omega=linspace(0,50,1024);

%% Covariance Matrix
C = zeros(length(tau));
for i = 1 : length(tau)
    for j=1:length(tau)
        C(i,j)= (mean_data(i)-data(i,:))*(mean_data(j)-data(j,:))';
    end
end
C =C/length(tau)/(length(tau)-1);
C_1=inv(C);
%%
figure(1)
subplot(1,3,1)
errorbar(tau,mean_data,(var_data))
title("imaginary time correlation function")
xlabel("\tau")
ylabel("c(\tau)")
subplot(1,3,2)
pcolor(C)
title("Covariance Matrix HeatMap")

%% Kernel
% integral of exp(-\tau\omega)*A(\omega)d\omega
integr = ones(length(omega),1)*2;
integr(1)=1;
integr(end)=1;
integr=integr*(omega(2)-omega(1))*0.5;

% compute the kernel
k= exp(- tau'.*omega);
% count the non singular eigenvalues
[V,S,U]=svd(k);
Ns=sum(sum((S>singular_value_tol)));

K=integr'.*k;
%% xhi_2
xhi_2 = @(A)(mean_data'-K*A)'* C_1*(mean_data'-K*A);

%% Information Entropy 
% flat distribution as prior
S = @(A) (omega(2)-omega(1))*sum(  A-1/(omega(end))-  (A>1e-5).*A.*log(A./(omega(end)))    );

%% Q 
% function to minimize for a givent alpha
Q = @(A,alpha) -alpha*S(A)+0.5*xhi_2(A);

get_A = @(u) exp(U(:,1:Ns)*u)/sum(exp(U(:,1:Ns)*u))/(omega(2)-omega(1));

%% Maxent Optimization

% For i in 1:nb alphas
logS=zeros(nb_alphas,1);
log_X = zeros(nb_alphas,1);
%options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxIter',500);
options = optimset('MaxIter',maxiter);
alphas=[alpha];
for i=1:nb_alphas 
    u=abs(normrnd(0,1,Ns,1));
    % Optimize alphaS-xhi_2/2 to find alpha
    u = fminsearch(@(x)Q(get_A(x),alpha),u,options);
    
    % save log(-S(A)) and log(\xhi_2)
    logS(i)=log(S(get_A(u)));
    log_X(i)=log(xhi_2(get_A(u)));
    % if i > 3 -> compute the curvature 
    figure(1)
    subplot(1,3,3)
    hold on
    plot(log_X,-logS,'x');
    title("L-Curve. Do not trust if it is not L-shaped");
    xlabel("log(\xhi)");
    ylabel("-log(S)");
    % Compute a new alpha
    alpha=alpha*step_alpha;
    alphas=[alphas alpha];
end
%% Find the minimum of curvature

f = zeros(nb_alphas-2,1);
[~,index]=sort(log_X);

for i = 1:length(f)
    f(i)= 2*logS(index(i+1))-logS(index(i))-logS(index(i+2));
end
[~,bestalpha]=max(abs(f));
%% Compute the definitive answer
alpha_fin=alphas(index(bestalpha+1));
u = fminsearch(@(x)Q(get_A(x),alpha_fin),u,options);

%% final plot
figure(2)
subplot(1,2,1)
plot(tau,mean_data)
hold on
plot(tau,K*get_A(u))
xlabel("\tau")
ylabel("C(\tau)")
legend("data","prediction")

subplot(1,2,2)
plot(omega,get_A(u))
xlabel("\omega")
ylabel("power Spectrum")

