function combMatrix = updateCombinationMatrix(combMatrix,numAgents,w,phi,gamma,A,attacker)
V = 0.1;    %%Forgetting factor

    for i = 1:numAgents
        invGammaSum = 0;  %%hold sum of inverse variance values in neighborhood of agent i
      if (i~= attacker)
        for j = 1:numAgents
            if (A(i,j) == 1 || (i == j))
                gamma(j,i) = (1-V)*gamma(j,i)+V*norm(phi(:,j)-w(:,i))^2; %%estimate variance of j relative to i
                invGammaSum = invGammaSum+ 1/gamma(j,i);
            end
        end
        for j = 1:numAgents         %%Once invGammaSum is calculated, update column i of combMatrix
            if (gamma(j,i) ~= 0)
                combMatrix(j,i) = (1/(gamma(j,i)))/invGammaSum;
            end
        end
    
      end
      
    
    end