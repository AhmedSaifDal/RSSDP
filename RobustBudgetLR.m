%% Description: Lagrangian Relaxation-Budgeted Uncertainty Set-G/M/1
clear
clc

% Inputs:
% I: The set of demand points
% J: The set of potential facility locations
% f: Setup cost (J)
% c: Access cost (I*J)
% xinom: Nominal demand (I)
% K: Number of breakpoints
tic
load('p170.mat');

f = 10*ones(J,1);
t = 10; % Waiting time penalty
CV = 0.5; % CV
R=(CV^2+1)/2;
xihat = Max_Dev;

K = 10;
Deviation = abs(bsxfun(@minus,SampleMat, xinom));
RelativeDeviation = bsxfun(@rdivide,Deviation, Max_Dev);
% Uncertainty Budget-Budgeted Uncertainty Set
RR = zeros(K,1);  
for q = 1:K
    d = sum(RelativeDeviation(q,:));
    RR(q,1) = d;
end

Gamma = prctile(R,90);
% Gamma = 0;
%% Piecewise Approximation-Generating Breaking Points
%g=@(q)q^2/(1-q);
%gderiv=@(q)(1/(1-q)^2)-1 ----> First Derivative of g
ghat=@(p)p^2/(1-p);%----> Linear Approximation
%T=@(q)q^2*(1+ghat(p)+p)-2*q*(ghat+p)+ghat(p)----> for finding q
%B=@(p)p^2*(1+gderiv(q))-p*(gderiv(q)+gderiv(q)*q-g(q)-e)+gderiv(q)*q-g(q)-e
p=0;
e=0.001;
bp=0; % Breaking Points
Tp=[]; % The set of Tangency Points
while p < 0.99
    r1=[1/(1-p)-e -2*p/(1-p)+2*e p^2/(1-p)-e];
    qq=roots(r1);
    for i=1:size(qq,1)
        if qq(i) > p
            q=qq(i);
        else
        end
    end
    Tp=[Tp q];
    r2=[1/(1-q)^2 e-2*q/(1-q)^2 q^2/(1-q)^2-e];
    pp=roots(r2);
    for n=1:size(pp,1)
        if pp(n) > q
            pnew=pp(n);
        else
        end
    end
    bp=[bp pnew];
    p=pnew;
end

K=length(bp);
ghatnew=zeros(1,K);
for k=1:K
ghatnew(:,k)=ghat(bp(k));
end

%% Subproblem [LSPj]:
upsilon = 500000*ones(I,1);
chi = Gamma/(I)*ones(I,1);
UB=Inf;
LB=-Inf;
Newcut=[];
RHNewcut=[];
iter = 0;
X=[];                                   % All the Generated Cut

while UB-LB > 0.01
iter = iter+1;
sol=zeros(J,I+3+K+1+I+1);               % Final Solution of All the Subproblems
beta=zeros(J,1);                        % A vector of all the Obj. values of Subproblems

% Order of variables: y_ij, mu_j, rho_j, theta_j,lambda_jk,\gamma_j,\eta_ij,w_j
% Equality Constraints:
Aeq1=[sparse(1,I), 0 , 1 , 0 , -bp, 0, sparse(1,I), 0];
beq1=0;

Aeq2=[sparse(1,I), 0 , 0 , 1 , -ghatnew, 0, sparse(1,I), 0];
beq2=0;

Aeq3=[sparse(1,I), 0 , 0 , 0 , ones(1,K), 0, sparse(1,I), 0];
beq3=1;

% w_j = Gamma*gamma_j+sum_i(eta_ij)
Aeq4=[sparse(1,I),0, 0, 0, sparse(1,K),-Gamma,-ones(1,I),1];
beq4=0;

Aeq = [Aeq1;Aeq2;Aeq3;Aeq4];
beq =[beq1;beq2;beq3;beq4];

% Inequality Constraints:
Aineq1=[diag(xihat),sparse(I,1),sparse(I,1),sparse(I,1),sparse(I,K),-ones(I,1),-speye(I,I),sparse(I,1)];
bineq1=zeros(I,1);

% SOC Constraint:
% Order of variables: y_ij, mu_j, rho_j, theta_j, lambda_jk,\gamma_j,\eta_ij,w_j
Q = sparse(I+3+K+1+I+1,I+3+K+1+I+1);
Q(1:I,1:I) = diag(xinom);
Q(I+1,I+2) = -1; 
L = zeros(I+3+K+1+I+1,1);
L(I+3+K+1+I+1,1)= 1 ;
modela.quadcon.Qc = Q;
modela.quadcon.q = L;
modela.quadcon.rhs = 0;

% SOS2
modela.sos.type = 2;
modela.sos.index = I+3+(1:K)';
modela.sos.wieght = (1:K)';

modela.A = [Aeq;Aineq1];
modela.rhs = [beq;bineq1];
modela.sense = [repmat('=',1,4),repmat('<',1,I)];
modela.vtype = [repmat('B',1,I),repmat('C',1,3+K+1+I+1)];
modela.lb = zeros(1,I+3+K+1+I+1);
modela.ub = [ones(1,I),Inf,1,Inf,ones(1,K),Inf(1,1+I+1)];
    
% Objective function
% For y_ij
A=c.*repmat(xinom',1,J);
B=(c.*repmat(xihat',1,J)).*chi;
C=repmat(upsilon,1,J);
D=reshape(A+B-C,1,I*J);

for j=1:J
    modela.obj=[D(1,(j-1)*I+1:I*j), f(j,1), t ,t*R ,zeros(1,K),0,zeros(1,I),0];
    modela.modelsence='Min';
    resultj = gurobi(modela);
    solj = resultj.x;
    sol(j,:) = solj';
    betaj = resultj.objval;
    beta(j,1) = betaj;
end

LB = max(LB,sum(beta)+sum(upsilon));
y = sol(:,1:I);
mu = sol(:,I+1);
rho = sol(:,I+2);
theta = sol(:,I+3);
lambda = sol(:,I+4:I+3+K);
H = [y,mu,rho,theta];
X=[X;H];

%% Master Problem-Kelly's Cutting Plane [DMP]:
% Order of variables: chi_i, beta_j, upsilon_i
% Generating Cuts:
BB = -(c'.*repmat(xihat,J,1)).*y;
Aineq2 = [BB ,speye(J,J), y];
AA = c'.*repmat(xinom,J,1);
P = [AA, f, t.*ones(J,1), t*R.*ones(J,1)].*H;
F = sum(P');
bineq2 = F';

Newcut=[Newcut;Aineq2];
RHNewcut=[RHNewcut;bineq2];

Aineq3 = [ones(1,I),sparse(1,J), sparse(1,I)];
bineq3 = Gamma;

Aineq4= [speye(I,I), sparse(I,J), sparse(I,I)];
bineq4 = ones(I,1);

modelb.obj = [zeros(1,I), -ones(1,J), -ones(1,I)];
modelb.modelsence='Min';
modelb.A =[Newcut;Aineq3;Aineq4];
modelb.rhs = [RHNewcut;bineq3;bineq4];
modelb.sense = repmat('<',1,J*iter+1+I);
modelb.vtype = repmat('C',1,I+J+I);
modelb.lb = [zeros(1,I),-Inf(1,J),zeros(1,I)];
modelb.ub = Inf(1,I+J+I);
result = gurobi(modelb);
solMaster = result.x;
valMaster = result.objval;

UB = -valMaster;
upsilon = solMaster(I+J+1:I+J+I);
fprintf('\niter=%2g, LB=%2.3f, UB=%2.3f \n\n', iter, LB, UB);
end

%% Dantzing-Wolfe Decomposition [MP]
% Order of variables:  alpha_j, omega, delta_i
% (chi_i):
Xtrapy=(X(:,1:I))';
SS = repmat(c.*repmat(xihat',1,J),1,iter-1).*Xtrapy;
Aineqdzw=[SS, -ones(I,1), -speye(I,I)];
bineqdzw=zeros(I,1);

% (upsilon_i):
Aeqdzw1=[Xtrapy, sparse(I,1), sparse(I,I)];
beqdzw1=ones(I,1);

% (beta_j):
Aeqdzw2=[repmat(speye(J,J),1,iter-1), sparse(J,1), sparse(J,I)];
beqdzw2=ones(J,1);

% Objective Function:
PP = repmat([AA, f, t.*ones(J,1), t*R.*ones(J,1)],iter-1,1).*X;
O = sum(PP');
modelc.obj = [O, Gamma, ones(1,I)] ;
modelc.modelsence='Min';
modelc.A =[Aineqdzw;Aeqdzw1;Aeqdzw2];
modelc.rhs = [bineqdzw;beqdzw1;beqdzw2];
modelc.sense = [repmat('<',1,I),repmat('=',1,I+J)];
modelc.vtype = repmat('C',1,size(X,1)+1+I);
modelc.lb = zeros(1,size(X,1)+1+I);
modelc.ub = [ones(1,size(X,1)),Inf(1,1+I)];
% modelc.ub = Inf(1,size(X,1)+1+I);
result = gurobi(modelc);
soldzw = result.x;
valdzw = result.objval;
alpha = soldzw(1:size(X,1));
omega = soldzw(size(X,1)+1);
delta = soldzw(size(X,1)+2:size(X,1)+1+I);


gap=(valdzw-LB)/LB;
timeElapsed = toc;

fprintf('valdzw=%2.3f, gap=%2.3f \n\n', valdzw, gap);

