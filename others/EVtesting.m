function EVtesting
load Graph.mat
I = eye(Graph.N);
m =5*ones(Graph.N,1);
EigenV = Graph.D;
EigenV(end,end) = 0;
Graph.Xtrue = 10+ Graph.Xtrue;
CovX    = zeros(Graph.N);
CovY    = zeros(Graph.N);
V       = rand(Graph.N);
Mvar    = 1;
S       = Graph.S(end:-1:1,:);

EstX = zeros(Graph.N,1);
Stepsize = 0.1/Graph.N;
theta = zeros(Graph.N,1);
loop = 1;
q = m+ Graph.S'*pinv(sqrt(EigenV))*sqrt(0.01)*randn(Graph.N,1); 
while(loop<50000)
    

[Xtrue,theta] = UpdateXtrue(m,q,theta);
Y = q+ sqrt(Mvar)* randn(Graph.N,1);

EstX = (1-Stepsize)*EstX + Stepsize*(Y);
CovX = CovX + ((Y-m)*(Y-m)');
CovY = CovY + (Y)*(Y)';

V = GF.EVtracking(V,CovY+I);

Ftrue = abs(S*Graph.Xtrue);
%Ftrue = SortingFreq(Ftrue);
F     = abs(V'*Graph.Xtrue);
%EstError(loop) = (1/Graph.N)*norm(Xtrue - EstX)^2;
Error(loop)  = (1/Graph.N)*norm(Ftrue-F)^2;
loop = loop+1;
end
figure(1)
clf
hold on
plot(10*log10(Error),'r');
%plot(10*log10(EstError),'g');
figure(2)
clf
hold on
stem(Ftrue(1:end-1),'r')
stem(F(1:end-1),'db')
legend('True Freq','Est Freq')
end

function [Xtrue,theta] = UpdateXtrue(m,q,theta)
theta = 0.9*theta + q;
Xtrue = m + theta;
end
function Ftemp = SortingFreq(F)
Ftemp = F;
n = size(F,1);
Ftemp(1:n-1) = F(2:n);
Ftemp(n)     = F(1);
end