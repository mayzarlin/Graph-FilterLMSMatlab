function w = EdgeWeight(Node1,Node2)
theta = -1;
w = exp(theta*(norm(abs(Node1) - (Node2))))^2;
end