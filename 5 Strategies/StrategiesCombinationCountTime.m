clear;
clc;

%% Proximity disk
Delta=3; % interaction distance 

network_no=2; %% or 2 or 3
[X_n,n,N]=load_network(network_no,Delta);
X_c=X_n;X_cs=X_n;X_cta=X_n;X_atc=X_n;
% n = dimension of each robot (n>1)
% N = total number of robots
% X = n x N vector containing the initial robot positions


%% Some numerical integration parameters
dt=0.003; % numerical steplength
Tf=3; % final time
numAgents = N-1;
t=0; 
iter0=2;
MSD_n = 0;MSD_c = 0;MSD_cs = 0;MSD_cta = 0;MSD_atc = 0;
a=-1;
b=1;

%Initialize regression sequence and weights for agents
finalno_iter=ceil(Tf/dt)+iter0;
u=zeros(finalno_iter,numAgents);
for k = 1:numAgents
u(:,k) = mvnrnd(0,.8,finalno_iter);
end
%Estimation Parameter
w0 = X_n(:,N);

step = 0.05;
iter=iter0;

v=0.1;
w_n = zeros(2,numAgents);w_c = w_n;w_cs = w_n;w_cta=w_n;w_atc=w_n;

fai_cs=zeros(2,numAgents);fai_atc=fai_cs;
d_cs=zeros(1,numAgents);d_cta=d_cs;d_atc=d_cs;
e_cs=zeros(1,numAgents);e_cta=e_cs;e_atc=e_cs;
ux_n=zeros(numAgents,2);ux_c=ux_n;ux_cs=ux_n;ux_cta=ux_n;ux_atc=ux_n;

% %% Save plot to Video
% 
% aviobj = VideoWriter('example3.avi');
% open(aviobj)

while (t<=Tf)
DX_n=zeros(n,N);DX_c=DX_n;DX_cs=DX_n;DX_cta=DX_n;DX_atc=DX_n;

sum_ci=0;
tcount_n=0;tcount_c=0;tcount_cs=0;tcount_cta=0;tcount_atc=0;
for i=1:numAgents
v=a+(b-a)*rand(1);
A=disk(X_cs,N,Delta);
    %% NonCoop
    tn=tic;
    ux_n = u(iter:-1:iter-1,i)';
    d_ni = ux_n*w0+v;
    w_n(:,i) = w_n(:,i) +step.*ux_n'*(d_ni - ux_n*w_n(:,i)); 
    tcount_n=tcount_n+toc(tn);
    %% CentralFusion
    tc=tic;
    ux_c = ux_n;
    d_ci = ux_c*w0+v; 
    sum_ci=sum_ci+ux_c'*(d_ci-ux_c*w_c(:,i));
    tcount_c=tcount_c+toc(tc);
    %% Consensus
    tcs=tic;
    fai_cs=zeros(2,numAgents);
    A0_cs=zeros(numAgents,numAgents);
    for j=1:numAgents
        N_neibor_cs=sum(A(:,j))+1;
        if A(j,i)==1 || i==j
          A0_cs(j,i)=1/N_neibor_cs;
        end
    end
    for j=1:numAgents
        fai_cs(:,i)=fai_cs(:,i)+A0_cs(j,i)*w_cs(:,j);
    end
    ux_cs(i,:) = ux_n;
    d_cs(i) = ux_cs(i,:)*w0+v;
    e_cs(i)=d_cs(i)-ux_cs(i,:)*w_cs(:,i);
    w_cs(:,i)=fai_cs(:,i)+step*ux_cs(i,:)'*e_cs(i);
    tcount_cs=tcount_cs+toc(tcs);
    %% CTA
    tcta=tic;
    fai_cta=zeros(2,numAgents);
    A0_cta=zeros(numAgents,numAgents);
    for j=1:numAgents
        N_neibor_cta=sum(A(:,j))+1;
        if A(j,i)==1  || i==j
          A0_cta(j,i)=1/N_neibor_cta;
        end
    end
    for j=1:numAgents
        fai_cta(:,i)=fai_cta(:,i)+A0_cta(j,i)*w_cta(:,j);
    end
    ux_cta(i,:) = ux_n;
    d_cta(i) = ux_cta(i,:)*w0+v;
    e_cta(i)=d_cta(i)-ux_cta(i,:)*fai_cta(:,i);
    w_cta(:,i)=fai_cta(:,i)+step*ux_cta(i,:)'*e_cta(i);
    tcount_cta=tcount_cta+toc(tcta);
    %% ATC
    tatc=tic;
    ux_atc(i,:) = ux_n;
    d_atc(i) = ux_atc(i,:)*w0+v;
    e_atc(i)=d_atc(i)-ux_atc(i,:)*w_atc(:,i);
    fai_atc(:,i)=w_atc(:,i)+step*ux_atc(i)*e_atc(i);
    A0_atc=zeros(numAgents,numAgents);
    for j=1:numAgents
        N_neibor_atc=sum(A(:,j))+1;
        if A(j,i)==1 || i==j 
          A0_atc(j,i)=1/N_neibor_atc;
        end
    end
    w_atc(:,i) = zeros(2,1);
    for j=1:numAgents
        w_atc(:,i)=w_atc(:,i)+A0_atc(j,i)*fai_atc(:,j);
    end
    tcount_atc=tcount_atc+toc(tatc);
end
%% CentralFusion
    tc=tic;
    w_c=w_c+step.*sum_ci/numAgents;
    tcount_c=tcount_c+toc(tc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Update time
t=t+dt;
for i=1:numAgents
    %% NonCoop
    tn=tic;
    DX_n(:,i)=DX_n(:,i)+(w_n(:,i)-X_n(:,i));
    X_n(:,i)=X_n(:,i)+dt.*DX_n(:,i);
    tcount_n=tcount_n+toc(tn);
    %% CentralFusion
    tc=tic;
    DX_c(:,i)=DX_c(:,i)+(w_c(:,i)-X_c(:,i));
    X_c(:,i)=X_c(:,i)+dt.*DX_c(:,i);
    tcount_c=tcount_c+toc(tc);
    %% Consensus
    tcs=tic;
    DX_cs(:,i)=DX_cs(:,i)+(w_cs(:,i)-X_cs(:,i));
    X_cs(:,i)=X_cs(:,i)+dt.*DX_cs(:,i);
    tcount_cs=tcount_cs+toc(tcs);
    %% CTA
    tcta=tic;
    DX_cta(:,i)=DX_cta(:,i)+(w_cta(:,i)-X_cta(:,i));
    X_cta(:,i)=X_cta(:,i)+dt.*DX_cta(:,i);
    tcount_cta=tcount_cta+toc(tcta);
    %% ATC
    tatc=tic;
    DX_atc(:,i)=DX_atc(:,i)+(w_atc(:,i)-X_atc(:,i));
    X_atc(:,i)=X_atc(:,i)+dt.*DX_atc(:,i);
    tcount_atc=tcount_atc+toc(tatc);
end

% %% Plot the solution every 10 iterations
% 
%   if (mod(iter,10)==0)
%     figure(1)
%     hold off
%     %% NonCoop
%     subplot(2,3,1)
%     hold off
%     for i=1:N
%       plot(X_n(1,i),X_n(2,i),'o');
%       xlabel('NonCoop')
%       hold on;
%       axis([-2*Delta,2*Delta,-Delta,2*Delta]);      
%     end
%     %% CentralFusion
%     subplot(2,3,2)
%     hold off
%     for i=1:N
%         plot(X_c(1,i),X_c(2,i),'o');
%         xlabel('CentralFusion')
%         hold on;
%         axis([-2*Delta,2*Delta,-Delta,2*Delta]);
%     end
%     %% Consensus
%     subplot(2,3,3)
%     hold off
%     for i=1:N
%         plot(X_cs(1,i),X_cs(2,i),'o');
%         xlabel('Consensus')
%         hold on;
%         axis([-2*Delta,2*Delta,-Delta,2*Delta]);
%     end
%     %% CTA
%     subplot(2,3,4)
%     hold off
%     for i=1:N
%         plot(X_cta(1,i),X_cta(2,i),'o');
%         xlabel('CTA')
%         hold on;
%         axis([-2*Delta,2*Delta,-Delta,2*Delta]);
%     end
%     %% ATC
%     subplot(2,3,5)
%     hold off
%     for i=1:N
%         plot(X_atc(1,i),X_atc(2,i),'o');
%         xlabel('ATC')
%         hold on;
%         axis([-2*Delta,2*Delta,-Delta,2*Delta]);
%     end
% 
% 
% 
%     fig = gcf;
%     exportgraphics(fig,'plot.png');
%     p1=imread('plot.png');
%     writeVideo(aviobj,p1);  
%     drawnow;
%   end

%   iter=iter+1;
end

% close(aviobj)
% for k = 1:numAgents
%     MSD_n = MSD_n + norm(w0-w_n(:,k))^2;
%     MSD_c = MSD_c + norm(w0-w_c(:,k))^2;
%     MSD_cs = MSD_cs + norm(w0-w_cs(:,k))^2;
%     MSD_cta = MSD_cta + norm(w0-w_cta(:,k))^2;
% end
% MSD_n = MSD_n/numAgents;
% MSD_c = MSD_c/numAgents;
% MSD_cs = MSD_cs/numAgents;
% MSD_cta = MSD_cta/numAgents;
% MSD_atc = MSD_atc + norm(w0-w_atc(:,k))^2;
% MSD_n
% MSD_c
% MSD_cs
% MSD_cta
% MSD_atc
tcount_n
tcount_c
tcount_cs
tcount_cta
tcount_atc
