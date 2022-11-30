function groupVel = ComputeGroupVelocity(combMatrix,w,vel,groupVel,numAgents)
v = 0.5;

for i = 1:numAgents
    theta(:,i) = zeros(2,1)
    for j = 1:numAgents
        theta(:,i) = theta(:,i) +combMatrix(j,i)*groupVel(:,j);
    end
end
for k = 1:numAgents
    groupVel(:,k) = theta(:,k)+v*(vel(:,k)-theta(:,k));
end
