close all
clear all
clc
set(0,'defaultTextInterpreter','latex');

%% Loading of the matrix

%M = csvread('../data/M_pre.csv');
M = csvread('../data/M_post_withGroup.csv'); %%%%% GROUP 
group = M(2:end,end); %%%%%%%%%%%%%% GROUP 
M = M(1:end,1:end-1); %%%%%%%%%%%%%% GROUP 

threshold = 20;
M = M(2:end,2:end);
Au = sparse(M); 
W = Au-diag(diag(Au)); % remove selfloops
W(W<threshold+1) = 0; % 0 for values <=20 and then weights for the others.
Au = double(W>threshold); % adjacency matrix with 0s and 1s
%W = W/max(max(W));

G = graph(Au);
[bins, binsize] = conncomp(G);

pos = find(sum(Au)~=0);
%group = group(pos); %%%%%%%%%%%%%%%% GROUP 
A = Au(pos,pos); % only the giant component
GC = graph(A); % its graph

figure(1)
plot(G, 'NodeColor','k','EdgeAlpha',0.1);

N = numnodes(GC);
L = numedges(GC);

disp(['Node: ' num2str(N) ' - Links: ' num2str(L)]) 

%% Comparison between strength of link and degree - no isolated vertex
Au=A;
s = full(sum(W,2)); % strength of link based on weighted matrix
s = s(s>0);
k = full(sum(Au,2)); % degree of a node based on adjacency matrix
k = k(k>0);

u = unique(k);
s_av = zeros(length(u),1);
for i = 1:length(u)
    s_av(i) = mean(s(k==u(i))); 
end

p = polyfit(log(u'),log(s_av'),1);
figure(2)
loglog(k,s,'.');
hold on
loglog(u,exp(p(2)+log(u)*p(1)),'r-');
hold off
grid
xlabel('k')
ylabel('s')
title('Strenght vs Degree')
% That is, the strength of a vertex is simply proportional to its degree, yielding an exponent ? = 1, 
% and the two quantities provide therefore the same information on the system.

disp(['<k>: ' num2str(mean(sum(Au,2)))])

%% Cluster coefficient

k = full(sum(A,2));
CC = zeros(1,N);
for i=find(k>1)'
    neighbors=find(A(i,:));
    Li=A(neighbors,neighbors);
    CC(i)= 2 * sum(sum(triu(Li))) / (k(i)*(k(i)-1));
end

disp(['<C> = ' num2str(mean(CC))]);

figure(3);
plot(k, CC, '.');
hold on; 
xlim([min(k),max(k)]);
grid on;
ylabel('C');
xlabel('k');
title('Clustering coefficient');

%% Degree distribution
figure(4)
subplot(2,1,1)
hist(k,100)
xlabel('k (degree)')
title('Degree and strength distribution') 
subplot(2,1,2)
hist(s,100)
xlabel('s (strength)')

% bimodal behaviour, difficult to fit with known distribution

%% EXTRACT THE DISTRIBUTION 
% Uncomment the following lines to see the weighted degree distribution
k = s;
u = unique(s);
pk = histc(k,u)'; % counts occurrences
pk = pk/sum(pk); % normalize to 1
Pk = cumsum(pk,'reverse'); % cumulative distribution
% log binning
klog = 10.^(0:0.1:ceil(log10(max(u))));
pklog = histc(k,klog)'; % counts occurrences
pklog = pklog/sum(pklog); % normalize to 1

figure(5)
subplot(2,2,1)
plot(u,pk,'.')
grid
xlabel('k')
ylabel('PDF')
title('linear PDF plot')
subplot(2,2,2)
loglog(u,pk,'.')
grid
xlabel('k')
ylabel('PDF')
title('logarithmic PDF plot')
subplot(2,2,3)
loglog(klog,pklog,'.')
grid
xlabel('k')
ylabel('PDF')
title('logarithmic PDF plot (log bins)')
subplot(2,2,4)
loglog(u,Pk,'.')
grid
xlabel('k')
ylabel('CCDF')
title('logarithmic CCDF plot')

%% PURE ML FITTING + ML FITTING WITH SATURATION
kmin = 600;
k2 = k(k>=kmin); % restrict range
ga = 1+1/mean(log(k2/kmin)); % estimate the exponent
disp(['gamma ML = ' num2str(ga)])

s1 = u.^(1-ga); % build the CCDF signal
s1 = s1/s1(150)*Pk(150);

for ks = 1:max(u)
    tmp = mean(log((k2+ks)/(kmin+ks)));
    ga2(ks) = 1+1/tmp;
    de(ks) = log(ga2(ks)-1)-log(kmin+ks)-ga2(ks)*tmp;
end
[~,ks] = max(de);
disp(['k_sat ML sat = ' num2str(ks)])
disp(['gamma ML sat = ' num2str(ga2(ks))])

s2 = ((u+ks)/(kmin+ks)).^(1-ga2(ks));
s2 = s2/s2(30)*Pk(30);

figure(6)
loglog(u,Pk,'.')
hold on
% ML fitting (we make sure that the plot follows the data)
loglog(u,s1);
% ML fitting with saturation
loglog(u,s2)
hold off
axis([xlim min(Pk/2) 2])
grid
xlabel('k')
ylabel('CCDF')
title('ML fittings')
legend('data','ML','ML with sat.')


%% other stuff
dist = distances(GC);
disp(['(threshold = ' num2str(threshold) ')']) 
disp(['Max distance (diamiter) = ' num2str(max(max(dist)))]) 
disp(['Average distance = ' num2str(mean(mean(dist)))]) 

k2 = k.^2;
disp(['<k^2>= ' num2str(mean(k2))]);
k3 = k.^3;
disp(['<k^3>= ' num2str(mean(k3))]);
k4 = k.^4;
disp(['<k^4>= ' num2str(mean(k4))]);

%% Assortativity
% only in the largest connect component (A inseed of Au)
% undirected network --> in-in/in-out/... assortativity is always the same 

d = sum(A,2); % degrees
k_tmp = (A*d)./d; % averages of neighbours

% extract averages for each value of k
u = unique(d);
k_nn = zeros(length(u),1);
for k = 1:length(u)
    k_nn(k) = mean(k_tmp(d==u(k)));
end



% linear fitting
p = polyfit(log(u),log(k_nn),1);
disp(['Assortativity factor = ' num2str(p(1))])

% cutoff
ks = (2*L)^0.5;
disp(['Structural cutoff = ' num2str(ks)])

% linear fitting with cutoff
p_cut = polyfit(log(u(u>514)),log(k_nn(u>514)),1);
disp(['Assortativity factor (with cutoff) = ' num2str(p_cut(1))])

% show
figure(7)
loglog(d,k_tmp,'g.');
hold on
loglog(u,exp(p(2)+log(u)*p(1)),'r-');
loglog(u,exp(p_cut(2)+log(u)*p_cut(1)),'b-');
loglog(u,k_nn,'k.');
hold off
grid
xlabel('k')
ylabel('$k_{nn}$')
title('Assortativity')



