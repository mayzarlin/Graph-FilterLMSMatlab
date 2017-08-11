% Node class contains properties and functions relating to each node in the Graph.
classdef RLMS < handle
    
    methods (Static)
 function node = GlobalAdapting(node,Graph)
            Stepsize = Graph.Stepsize/Graph.N;
        %% Adapting
         node.GlobalX =  (1-Stepsize)*node.GlobalX + Stepsize*Graph.Y; 
        end
 function node = GlobalReg(node,Graph)
      I = eye(Graph.N); L = Graph.Lsym;
      beta = 2*((node.GlobalX'*L*node.GlobalX)-Graph.Xtrue'*L*Graph.Xtrue)/norm(L*node.GlobalX)^2;
      beta = max(beta,0);
      beta = min(beta,1);
       
%% Filtering
       node.GlobalX = (I-(beta*L))*node.GlobalX;
 end
 function node = LocalReg(node,Graph)
beta = zeros(Graph.N,1);
I = eye(Graph.N);
                %% Local Beta
                for i = 1:Graph.N
                    Local_Reg           = 2*((node.X(:,i)'*Graph.Local(i).L* node.X(:,i))....
                                                -(Graph.Xtrue'*Graph.Local(i).L* Graph.Xtrue));
                    Local_Regnorm       =   norm( Graph.Local(i).L*node.X(:,i))^2;
                    beta_local          = Local_Reg/Local_Regnorm;
                    beta_local          = min(beta_local,1);
                    beta(i)  = max(beta_local,0);
                end
                %% Filtering
                for i = 1:Graph.N
                    node.X(:,i)            = (I - beta(i)*Graph.Local(i).L)*node.X(:,i);
                end
end
end
end


