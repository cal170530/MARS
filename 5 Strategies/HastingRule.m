function[fai]=HastingRule(gamma,w,v,i,numAgents,A)
sum_gamma=0;
fai=zeros(2,numAgents);
A0=zeros(numAgents,numAgents);
for j=1:numAgents
    if A(j,i)==1
      gamma(j,i)=(1-v)*gamma(j,i)+v*norm(fai(:,i)-w(:,i))^2;
      if gamma(i,j)~=0
       sum_gamma=sum_gamma+gamma(j,i)^(-1);
      end
    end
end
for j=1:numAgents
    if A(j,i)==1 && sum_gamma~=0
      A0(j,i)=1/sum_gamma/gamma(j,i);
    end
    fai(:,i)=fai(:,i)+A0(j,i)*w(:,j);
end
