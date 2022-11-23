%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

figure
%% Proximity disk
Delta=3; % interaction distance

network_no=4; %% or 2 or 3
[X,n,N]=load_network(network_no,Delta);
% n = dimension of each robot (n>1)
% N = total number of robots
% X = n x N vector containing the initial robot positions


%% Some numerical integration parameters
dt=0.003; % numerical steplength
Tf=5; % final time
numAgents = N-1;
t=0; 
iter0=2;
MSD = 0;

% 1.Jingyi comment:the changes below is just to simplify codes and increase adaptability
w = zeros(2,numAgents);
%Initialize regression sequence and weights for agents
finalno_iter=ceil(Tf/dt)+iter0;
for k = 1:numAgents
u(:,k) = mvnrnd(0,.8,finalno_iter);
end
%Estimation Parameter
w0 = X(:,N);

step = 0.05;
iter=iter0;
while (t<=Tf)
% 2.Jingyi comments:note that the w in NonCoop case is calculated by each
% single nodes. And each i node can update it's w(:,i) and X(:,i) seperately

DX=zeros(n,N); %% Here is where we store the derivatives 

for i=1:numAgents
ux = u(iter:-1:iter-1,i)';
d_i = ux*w0+normrnd(0,.15);
w(:,i) = w(:,i) +step.*ux'*(d_i - ux*w(:,i)); 
%%%%%%%%%%%%%%%%%%%%%  
DX(:,i)=DX(:,i)+(w(:,i)-X(:,i));
% Update the states using an Euler approximation
X(:,i)=X(:,i)+dt.*DX(:,i);
end


%% Update time
  t=t+dt;


%% Plot the solution every 10 iterations
% 3. Jingyi comments:here I add a state code 0 to plot without any edge.
% With state == 1, there will be edges connected every two nodes.
  if (mod(iter,10)==0)
    plotsol(X,N,Delta,0);
  end

  iter=iter+1;
end

for k = 1:numAgents
    MSD = MSD + norm(w0-w(:,k))^2;
end
MSD = MSD/numAgents;
MSD
%4. Jingyi comments:MSD_ncop is approximate to step*M*sigma^2/2. 
% Note that MSD_ncop is around N times as MSD_cent


