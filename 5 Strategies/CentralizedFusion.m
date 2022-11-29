%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

figure
%% Proximity disk
Delta=3; % interaction distance

network_no=2; %% or 2 or 3
[X,n,N]=load_network(network_no,Delta);
% n = dimension of each robot (n>1)
% N = total number of robots
% X = n x N vector containing the initial robot positions


%% Some numerical integration parameters
dt=0.003; % numerical steplength
Tf=3; % final time
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
% 2.Jingyi comments:note that the w in CentralizedFussion case is without subscript k, meaning that
% it's the estimator calculated by the whole network
sum_i=0;
for i=1:numAgents
    ux = u(iter:-1:iter-1,i)';
    d_i = ux*w0+normrnd(0,.15); 
    sum_i=sum_i+ux'*(d_i-ux*w(:,i));
end
w=w+step.*sum_i/numAgents;
%% Update the states using an Euler approximation
DX=zeros(n,N); %% Here is where we store the derivatives
for i=1:numAgents
    DX(:,i)=DX(:,i)+(w(:,i)-X(:,i));
    X(:,i)=X(:,i)+dt.*DX(:,i);
end

%% Update time
t=t+dt;


%% Plot the solution every 10 iterations
% 3. Jingyi comments:here I add a state code 1 to plot with edges connected every two nodes.
% Otherwise,it will plot without edges
if (mod(iter,10)==0)
    plotsol(X,N,Delta,1);
end

iter=iter+1;
end

for k = 1:numAgents
    MSD = MSD + norm(w0-w(:,k))^2;
end
MSD = MSD/numAgents;
MSD
% 4. Jingyi comments: MSD_cent is approcimate to step*M*sigma^2/2/numAgents
% Note that MSD_ncop is around N times as MSD_cent

