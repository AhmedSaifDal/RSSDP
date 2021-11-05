%% Description: ROSSDP-Budgeted Uncertainty Set-G/M/1
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
load('p2.mat');

% f=60*ones(J,1);
t=10; % Waiting time penalty
CV=0.5; % CV
R=(CV^2+1)/2;
xihat=Max_Dev;
% Gamma = 0;
params.MIPGap = 0.001;
params.TIME_LIMIT = 10000;
% params.NonConvex=2;
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

%% Problem
% Order of variables: y_ij, mu_j, rho_j, theta_j,lambda_jk,\theta',\phi'_i,\gamma_j,\eta_ij,w_j
% Equality Constraints:
Aeq1=[repmat(speye(I,I),1,J),sparse(I,J), sparse(I,J), sparse(I,J), sparse(I,J*K), sparse(I,1), sparse(I,I),sparse(I,J),sparse(I,I*J),sparse(I,J)];
beq1=ones(I,1);

% bp: Breakpoints (1*K vector)
Aeq2=[sparse(J,I*J), sparse(J,J), speye(J,J), sparse(J,J), -kron(speye(J,J),bp), sparse(J,1), sparse(J,I),sparse(J,J),sparse(J,I*J),sparse(J,J)];
beq2=zeros(J,1);

Aeq3=[sparse(J,I*J), sparse(J,J), sparse(J,J), speye(J,J), -kron(speye(J),ghatnew), sparse(J,1), sparse(J,I),sparse(J,J),sparse(J,I*J),sparse(J,J)];
beq3=zeros(J,1);

Aeq4=[sparse(J,I*J),sparse(J,J), sparse(J,J), sparse(J,J), kron(speye(J),ones(1,K)), sparse(J,1), sparse(J,I),sparse(J,J),sparse(J,I*J),sparse(J,J)];
beq4=ones(J,1);

Aeq5=[sparse(J,I*J),sparse(J,J), sparse(J,J), sparse(J,J), sparse(J,J*K), sparse(J,1), sparse(J,I),-Gamma*speye(J,J),-kron(eye(J,J),ones(1,I)),speye(J,J)];
beq5=zeros(J,1);

Aeq = [Aeq1;Aeq2;Aeq3;Aeq4;Aeq5];
beq =[beq1;beq2;beq3;beq4;beq5];

% Inequality Constraints:
xx=c.*repmat(xihat',1,J);
X=[];
for i=1:size(xx,2)
    ss=diag(xx(:,i));
    X=[X ss];
end
Aineq1=[X,sparse(I,J),sparse(I,J),sparse(I,J),sparse(I,J*K),-ones(I,1),-speye(I,I),sparse(I,J),sparse(I,I*J),sparse(I,J)];
bineq1=zeros(I,1);

Aineq2=[kron(speye(J),diag(xihat')),sparse(I*J,J),sparse(I*J,J),sparse(I*J,J),sparse(I*J,J*K),sparse(I*J,1),sparse(I*J,I),kron(speye(J,J),-ones(I,1)),-speye(I*J,I*J),sparse(I*J,J)];
bineq2=zeros(I*J,1);

Aineq = [Aineq1;Aineq2];
bineq =[bineq1;bineq2];

model.A = [Aeq;Aineq];
model.rhs = [beq;bineq];

% SOC Constraint:
% Order of variables: y_ij, mu_j, rho_j, theta_j, lambda_jk,\theta',\phi'_i,\gamma_j,\eta_ij,w_j
for j=1:J
L = zeros(I*J+3*J+J*K+1+I+J+I*J+J,1);
L(I*J+3*J+J*K+1+I+J+I*J+j,1)=1;
Q = sparse(I*J+3*J+J*K+1+I+J+I*J+J,I*J+3*J+J*K+1+I+J+I*J+J);
Q((j-1)*I+1:j*I,(j-1)*I+1:j*I)= diag(xinom);
Q(I*J+j,I*J+J+j)= -1; 
model.quadcon(j).Qc = Q;
model.quadcon(j).q = L;
model.quadcon(j).rhs = 0;
model.quadcon(j).name = sprintf('cone%d',j);
end

% SOS2
for i = 1:J
    model.sos(i).type = 2;
    model.sos(i).index = I*J+3*J+K*(i-1)+[1:K]';
    model.sos(i).wieght = [1:K]';
end

% Objective function
% For y_ij
D=reshape(c.*repmat(xinom',1,J),1,I*J);
model.obj=[D,f',t.*ones(1,J),t*R.*ones(1,J),zeros(1,J*K),Gamma,ones(1,I),zeros(1,J),zeros(1,I*J),zeros(1,J)];
model.modelsence='Min';
model.sense = [repmat('=',1,I+3*J),repmat('<',1,I+I*J+J)];

model.vtype = [repmat('B',1,I*J),repmat('C',1,3*J+J*K+1+I+J+I*J+J)];
model.lb = zeros(1,I*J+3*J+J*K+1+I+J+I*J+J);
model.ub = [ones(1,I*J),Inf(1,J),ones(1,J),Inf(1,J),ones(1,J*K),Inf(1,1+I+J+I*J+J)];

gurobi_write(model, 'qcp.lp');
result = gurobi(model, params);

sol = result.x;
val = result.objval;
y=reshape(sol(1:I*J),I,J);
mu=sol(I*J+1:I*J+J);
rho=sol(I*J+J+1:I*J+2*J);
lambda=reshape(sol(I*J+3*J+1:I*J+3*J+J*K),J,K);
% save ROBud270 y mu rho lambda

timeElapsed = toc;

val
timeElapsed
