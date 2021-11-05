%% Description: ROSSDP-Ball Uncertainty Set-G/M/1
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
load('p1.mat');

t=10; % Waiting time penalty
CV=0; % CV
R=(CV^2+1)/2;
xihat=Max_Dev;
% r = rhoBall;
r = rho;
params.MIPGap = 0.001;
params.TIME_LIMIT = 10000;
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
    r1=[e-1/(1-p) 2*p/(1-p)-2*e e-p^2/(1-p)];
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
% Order of variables: y_ij, mu_j, rho_j, theta_j, lambda_jk,u ,u'_j
% Equality Constraints:
Aeq1=[repmat(speye(I,I),1,J),sparse(I,J), sparse(I,J), sparse(I,J), sparse(I,J*K), sparse(I,1), sparse(I,J)];
beq1=ones(I,1);

% bp: Breakpoints (1*K vector)
Aeq2=[sparse(J,I*J), sparse(J,J), speye(J,J), sparse(J,J), -kron(speye(J,J),bp),sparse(J,1),sparse(J,J)];
beq2=zeros(J,1);

Aeq3=[sparse(J,I*J), sparse(J,J), sparse(J,J), speye(J,J), -kron(speye(J),ghatnew),sparse(J,1),sparse(J,J)];
beq3=zeros(J,1);

Aeq4=[sparse(J,I*J),sparse(J,J), sparse(J,J), sparse(J,J), kron(speye(J),ones(1,K)),sparse(J,1),sparse(J,J)];
beq4=ones(J,1);

Aeq = [Aeq1;Aeq2;Aeq3;Aeq4];
beq =[beq1;beq2;beq3;beq4];

model.A = Aeq;
model.rhs = beq;

% SOC Constraints:
%(1):
Q1=sparse(I*J+3*J+J*K+1+J,I*J+3*J+J*K+1+J);
A=(reshape(c,1,I*J)).^2;
Q1(1:I*J,1:I*J)=diag(A);
Q1(I*J+3*J+J*K+1,I*J+3*J+J*K+1)=-1; 
model.quadcon(1).Qc= Q1;
model.quadcon(1).q = zeros(I*J+3*J+J*K+1+J,1);
model.quadcon(1).rhs = 0;
model.quadcon(1).name = sprintf('cone%d',1);

%(2):
for j=1:J
L = zeros(I*J+3*J+J*K+1+J,1);
L(I*J+3*J+J*K+1+j,1)= r;
Q2=sparse(I*J+3*J+J*K+1+J,I*J+3*J+J*K+1+J);
Q2((j-1)*I+1:j*I,(j-1)*I+1:j*I)=diag(xinom);
Q2(I*J+j,I*J+J+j)=-1; 
model.quadcon(1+j).Qc= Q2;
model.quadcon(1+j).q = L;
model.quadcon(1+j).q = zeros(I*J+3*J+J*K+1+J,1);
model.quadcon(1+j).rhs = 0;
model.quadcon(1+j).name = sprintf('cone%d',1+j);
end

%(3):
for j=1:J
Q3=sparse(I*J+3*J+J*K+1+J,I*J+3*J+J*K+1+J);
Q3((j-1)*I+1:j*I,(j-1)*I+1:j*I)=diag(ones(1,I));
Q3(I*J+3*J+J*K+1+j,I*J+3*J+J*K+1+j)=-1; 
model.quadcon(1+J+j).Qc= Q3;
model.quadcon(1+J+j).q = zeros(I*J+3*J+J*K+1+J,1);
model.quadcon(1+J+j).rhs = 0;
model.quadcon(1+J+j).name = sprintf('cone%d',1+J+j);
end


% SOS2
for i = 1:J
    model.sos(i).type = 2;
    model.sos(i).index = I*J+3*J+K*(i-1)+[1:K]';
    model.sos(i).wieght = [1:K]';
end

% Objective function
% Order of variables: y_ij, mu_j, rho_j, theta_j, lambda_jk,u,u'_j
% For y_ij
D=reshape(c.*repmat(xinom',1,J),1,I*J);
model.obj=[D,f',t.*ones(1,J),t*R.*ones(1,J),zeros(1,J*K),r,zeros(1,J)];
model.modelsence='Min';
model.sense = repmat('=',1,I+3*J);

model.vtype = [repmat('B',1,I*J),repmat('C',1,3*J+J*K+1+J)];
model.lb = zeros(1,I*J+3*J+J*K+1+J);
model.ub = [ones(1,I*J),Inf(1,J),ones(1,J),Inf(1,J),ones(1,J*K),Inf(1,1+J)];

gurobi_write(model, 'qcp.lp');
result = gurobi(model, params);

sol = result.x;
val = result.objval;
y=reshape(sol(1:I*J),I,J);
mu=sol(I*J+1:I*J+J);
rho=sol(I*J+J+1:I*J+2*J);
lambda=reshape(sol(I*J+3*J+1:I*J+3*J+J*K),J,K);
% save ROBall270 y mu rho lambda

timeElapsed = toc;

val
timeElapsed

