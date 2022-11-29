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

%Initialize regression sequence and weights for agents
finalno_iter=ceil(Tf/dt)+iter0; 
for k = 1:numAgents
    u(:,k) = mvnrnd(0,.8,finalno_iter);
end
%Estimation Parameter
w0 = X(:,N);


step = 0.05;
iter=iter0;

v=0.1;
w = zeros(2,numAgents);
fai=zeros(2,numAgents);
d=zeros(1,numAgents);
e=zeros(1,numAgents);
ux=zeros(numAgents,2);
gamma=zeros(numAgents,numAgents);
% A=zeros(numAgents,numAgents);
tic
while (t<=Tf)
% 2.Jingyi comments:implement Hasting rule for left-stochastic combination
% matrices(eq. 33)

for i=1:numAgents
    ux(i,:) = u(iter:1:iter+1,i)';
    d(i) = ux(i,:)*w0+normrnd(0,.15);
    e(i)=d(i)-ux(i,:)*w(:,i);
    fai(:,i)=w(:,i)+step*ux(i)*e(i);
% end

%%     Hasting's Rule
% for i=1:numAgents
    A=disk(X,N,Delta);
    sum_gamma=0;
    A0=zeros(numAgents,numAgents);
    for j=1:numAgents
        if A(j,i)==1
          gamma(j,i)=(1-v)*gamma(j,i)+v*norm(fai(:,j)-w(:,i))^2;
          if gamma(i,j)~=0
            sum_gamma=sum_gamma+1/gamma(j,i);
          end
        end
    end
% end


    w(:,i) = zeros(2,1);
%     A=disk(X,N,Delta);
    for j=1:numAgents
        if A(j,i)==1 && sum_gamma~=0
          A0(j,i)=1/sum_gamma/gamma(j,i);
          
        end
        w(:,i)=w(:,i)+A0(j,i)*fai(:,j);
    end
    %%
end
    


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
toc
for k = 1:numAgents
    MSD = MSD + norm(w0-w(:,k))^2;
end
MSD = MSD/numAgents;
MSD
% 4. Jingyi comments: MSD_cent is approcimate to step*M*sigma^2/2/numAgents
% Note that MSD_ncop is around N times as MSD_cent

