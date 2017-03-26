function OmegaVal=LMI_Aut16b_th2(A,B,beta,D0,Q,C,K,rho,BC,Gamma,tauM,Delta,epsilon,alpha)
% This MATLAB program checks the feasibility of LMIs from Theorem 2 of the paper 
% A. Selivanov and E. Fridman, "Distributed event-triggered control of diffusion semilinear PDEs," Automatica, vol. 68, pp. 344–351, 2016.

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and SeDuMi solver (http://sedumi.ie.lehigh.edu/)

% Input: 
% A,B,beta      - the parameters of the system (1); 
% D0            =diag{d_1^0,...,d_n^0}, 0<d_i^0<=d_i(x), where d_i(x) is from (1); 
% Q             - the nonlinearity estimate in the form of (2); 
% C             - the measurements matrix from (6); 
% K             - the controller gain from (10); 
% rho           - diagonal matrix of the quantizers' densities; 
% BC            - the type of boundary conditions: 1 - Dirichlet (3), 2 - Neumann (4), 3 - mixed (5); 
% Gamma         - the matrix from (5); 
% tauM          - h+MAD, the overall delay; 
% Delta         - the maximum size of the subdomains; 
% epsilon       - the event-triggering parameter from (8); 
% alpha         - the decay rate; 

% Output: 
% OmegaVal      - the event-triggering matrix from (8). If OmegaVal is empty, the LMIs are not feasible. 

[m,n]=size(C); 

%% Decision variables 
P1=sdpvar(n); 
P2=diag(sdpvar(1,n)); 
P3=diag(sdpvar(1,n)); 
R=sdpvar(n); 
S=sdpvar(n); 
Omega=sdpvar(m); 
Lambdaq=diag(sdpvar(1,m)); 
Lambdatheta=diag(sdpvar(1,n)); 
Lambdakappa=diag(sdpvar(1,n)); 
G=sdpvar(n,n,'f');
lambdaphi=sdpvar; 

%% Notations
Deltaq=diag((1-rho)./(1+rho)); 

%% The LMI for Psi
Psi=blkvar; 
Psi(1,1)=S-exp(-alpha*tauM)*R+P2*A+A'*P2-P2*B*K*C+alpha*P1-(P2*B*K*C)'+...
    lambdaphi*Q+Lambdatheta+C'*Lambdaq*Deltaq^2*C+epsilon*C'*Omega*C; 
Psi(1,2)=P1-P2+A'*P3-(P3*B*K*C)'; 
if BC==2 || ~isdiag(beta)
	Psi(1,3)=-P2*beta; 
end
Psi(1,4)=exp(-alpha*tauM)*G'; 
Psi(1,5)=exp(-alpha*tauM)*(R-G')-Lambdatheta; 
Psi(1,6)=P2; 
Psi(1,7)=-P2*B*K*C+C'*Lambdaq*Deltaq^2*C+epsilon*C'*Omega*C; 
Psi(1,8)=-P2*B*K*C+C'*Lambdaq*Deltaq^2*C+epsilon*C'*Omega*C;
Psi(1,9)=-P2*B*K; 
Psi(1,10)=-P2*B*K+C'*Lambdaq*Deltaq^2; 
Psi(2,2)=tauM^2*R-2*P3; 
Psi(2,3)=-P3*beta; 
Psi(2,6)=P3; 
Psi(2,7)=-P3*B*K*C; 
Psi(2,8)=-P3*B*K*C; 
Psi(2,9)=-P3*B*K; 
Psi(2,10)=-P3*B*K; 
Psi(3,3)=D0*(alpha*P3-2*P2)+Delta^2/pi^2*Lambdakappa; 
Psi(4,4)=-exp(-alpha*tauM)*(S+R); 
Psi(4,5)=exp(-alpha*tauM)*(R-G); 
Psi(5,5)=-2*exp(-alpha*tauM)*R+exp(-alpha*tauM)*(G+G')+Lambdatheta; 
Psi(6,6)=-lambdaphi*eye(n); 
Psi(7,7)=-Lambdatheta+C'*Lambdaq*Deltaq^2*C+epsilon*C'*Omega*C; 
Psi(7,8)=C'*Lambdaq*Deltaq^2*C+epsilon*C'*Omega*C; 
Psi(7,10)=C'*Lambdaq*Deltaq^2; 
Psi(8,8)=-Lambdakappa+C'*Lambdaq*Deltaq^2*C+epsilon*C'*Omega*C; 
Psi(8,10)=C'*Lambdaq*Deltaq^2; 
Psi(9,9)=-Lambdaq; 
Psi(10,10)=Lambdaq*Deltaq^2-Omega; 
Psi=sdpvar(Psi); 

%% Park's condition
Park=[R G; G' R];

%% Solution of LMIs
LMIs=[P1>=0, P3>=0, R>=0, S>=0, Omega>=0, Lambdaq>=0, Lambdatheta>=0, Lambdakappa>=0, lambdaphi>=0, Psi<=0, Park>=0]; 
if BC==3 % The condition from (iii)
    LMIs=[LMIs, 2*(alpha*P3-2*P2)*D0*Gamma+P2*beta+beta'*P2<=0]; 
end
options = sdpsettings('solver','sedumi','verbose',0);
sol=optimize(LMIs,[],options); 

OmegaVal=[]; 
if sol.problem == 0
    primal=check(LMIs); 
    if min(primal)>=0 && min(primal(1:4))>0
        OmegaVal=value(Omega); 
    end
else
    yalmiperror(sol.problem) 
end
