function OmegaVal=LMI_Aut16b_th1(A,B,beta,D0,Q,C,K,rho,BC,Gamma,tauM,Delta,epsilon,delta,delta1)
% This MATLAB program checks the feasibility of LMIs from Theorem 1 of the paper 
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
% Gamma         - the matrix from (5). 
% tauM          - h+MAD, the overall delay; 
% Delta         - the maximum size of the subdomains; 
% epsilon       - the event-triggering parameter from (8); 
% delta, delta1 - determine the decay rate via (19); 

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
G=sdpvar(n,n,'f');
lambdaphi=sdpvar; 

%% Notations
Deltaq=diag((1-rho)./(1+rho)); 

%% The LMI for Xi
Xi=blkvar; 
Xi(1,1)=S-exp(-delta*tauM)*R+P2*A+A'*P2+lambdaphi*Q+delta*P1; 
Xi(1,2)=P1-P2+A'*P3; 
if BC==2 || ~isdiag(beta)
	Xi(1,3)=-P2*beta; 
end
Xi(1,4)=exp(-delta*tauM)*G'; 
Xi(1,5)=exp(-delta*tauM)*(R-G')-P2*B*K*C; 
Xi(1,6)=P2; 
Xi(1,7)=-P2*B*K*C; 
Xi(1,8)=-P2*B*K; 
Xi(1,9)=-P2*B*K;
Xi(2,2)=tauM^2*R-2*P3; 
Xi(2,3)=-P3*beta; 
Xi(2,5)=-P3*B*K*C; 
Xi(2,6)=P3; 
Xi(2,7)=-P3*B*K*C; 
Xi(2,8)=-P3*B*K; 
Xi(2,9)=-P3*B*K; 
Xi(3,3)=D0*(delta*P3-2*P2); 
Xi(4,4)=-exp(-delta*tauM)*(S+R); 
Xi(4,5)=exp(-delta*tauM)*(R-G); 
Xi(5,5)=-2*exp(-delta*tauM)*R+exp(-delta*tauM)*(G+G')+C'*Lambdaq*Deltaq^2*C+epsilon*C'*Omega*C-delta1*P1; 
Xi(5,7)=C'*Lambdaq*Deltaq^2*C+epsilon*C'*Omega*C; 
Xi(5,9)=C'*Lambdaq*Deltaq^2; 
Xi(6,6)=-lambdaphi*eye(n);
Xi(7,7)=C'*Lambdaq*Deltaq^2*C+epsilon*C'*Omega*C-delta1*P3*D0*pi^2/Delta^2; 
Xi(7,9)=C'*Lambdaq*Deltaq^2; 
Xi(8,8)=-Lambdaq; 
Xi(9,9)=Lambdaq*Deltaq^2-Omega; 
Xi=sdpvar(Xi); 

%% Park's condition
Park=[R G; G' R];

%% Solution of LMIs
LMIs=[P1>=0, P3>=0, R>=0, S>=0, Omega>=0, Lambdaq>=0, lambdaphi>=0, Xi<=0, Park>=0]; 
if BC==3 % The condition from (iii)
    LMIs=[LMIs, 2*(delta*P3-2*P2)*D0*Gamma+P2*beta+beta'*P2<=0]; 
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
