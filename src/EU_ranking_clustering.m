close all
clear all
clc
set(0,'defaultTextInterpreter','latex');

%% Loading of the matrix
% M = csvread('../data/M_pre.csv');
M = csvread('../data/M_post_withGroup.csv'); %%%%% GROUP 
group = M(2:end,end); %%%%%%%%%%%%%% GROUP 
M = M(1:end,1:end-1); %%%%%%%%%%%%%% GROUP 
names = readtable('../data/names.xlsx','ReadVariableNames', false);

threshold = 20;
M = M(2:end,2:end);
Au = sparse(M); 
W = Au-diag(diag(Au)); % remove selfloops

W(W<threshold+1) = 0; % 0 for values <=20 and then weights for the others.
Au = double(W>threshold); % adjacency matrix with 0s and 1s
% W = W/max(max(W));
G = graph(Au);
[bins, binsize] = conncomp(G);

pos = find(sum(Au)~=0);

group = group(pos); %%%%%%%%%%%%%%%% GROUP 

A = Au(pos,pos); % only the giant component
GC = graph(A);

index = 1:854;
index = index';
G.Nodes.Name = cellstr(num2str(index));
GC.Nodes.Name = cellstr(num2str(index(pos)));

%% Graph Plot
figure(1)
subplot(2,1,1)
plot(G, 'NodeColor','k','EdgeAlpha',0.1);
title(['Graph "as it is" with threshold=' num2str(threshold)])
subplot(2,1,2)
H = plot(GC, 'NodeColor','k','EdgeAlpha',0.1);
title('Giant component with label')
%H.NodeLabel = index(pos);
names_of_MEPs = table2array(names(pos,:));
H.NodeLabel = names_of_MEPs;


L = numedges(GC);
nnode = numnodes(GC);
disp(['Node: ' num2str(nnode) ' - Edge: ' num2str(L)])

%%  Pre-computations
d = full(degree(GC));
%A = W(pos,pos);
N = size(A,1);
D = diag(d.^-0.5); % degree matrix ^-1/2
L1 = spdiags(ones(N,1),0,N,N) - (D*A*D); % normalized laplacian
% u=eigenvector, lambda is a diagonal matrix
[u,lambda] = eigs(L1, N, 'smallestabs');  

%% %%%%%%%%%%%%%%%%%%%%%%%%% PageRank %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('PageRank')
% pre-computations and parameters
d = sum(A);
M = sparse(A*diag(1./d));  
c = 0.85;
q = repmat(1/N,N,1); 

% linear system equation
B = (sparse(eye(N))-c*M);
y = (1-c)*q;
disp('--- Linear system')
tic
pl = B\y;
toc

% power iteration method
p = q;
k = 40;  % max_iterations
%k = size(A,1);
eig_2 = eigs(M,size(M,1));
e = zeros(1,k);
r = zeros(1,k);
disp('--- Power iteration')
tic
for i=1:k
    r(1,i) = norm(p-pl);
    p = c*M*p + (1-c)*q;
    e(1,i) = (c*abs(eig_2(2)))^i;
end
toc

%% Visualization PageRank
figure(2)
semilogy(1:k, r)
hold on
semilogy(1:k, e/e(20)*r(20))
hold off
grid
xlabel('k')
ylabel('$||p-p_l||$')
title('PageRank convergence')
legend('power iteration','second eigenvalue')

figure(3)
plot(eig_2,'x')
hold on
plot(exp(2i*pi*(0:0.001:1)))
hold off
grid
title('PageRank eigenvalues')

%% extract MEPs with higher values
clc

k = 10;
p_sort = sort(p, 'descend');
p_sort = p_sort(1:k);
for i=1:k
    max_p_pos_s = find(p==p_sort(i));
    disp(names(max_p_pos_s,1));
    disp([' -- ranking value: ' num2str(p_sort(i))]);
end

%% HITS - authorities

M = A'*A;
disp('Computing HITS - eigenvalue extraction')
tic
[pp,ee] = eigs(M,2);
toc
p_hits = -pp(:,1)/norm(pp(:,1));


%% extract MEPs with higher values

k = 8;
p_sort = sort(p_hits, 'ascend');
p_sort = p_sort(1:k);
for i=1:k
    max_p_pos_s = find(p_hits==p_sort(i));
    disp(names(max_p_pos_s,1));
    disp([' -- ranking value: ' num2str(-p_sort(i),6)]);
end


%% comparison 

figure(4)
plot([p_hits/sum(p_hits),-p/sum(p)])
grid
legend('HITS','PageRank')
title('PageRank vs HITS ')

figure(5)
plot(p_hits/sum(p_hits),p/sum(p),'x')
grid
xlabel('HITS score')
ylabel('Pagerank score')
title('PageRank vs HITS ')



%% %%%%%%%%%%%%%%%% Spectral Approch %%%%%%%%%%%%%%%%%%%%%%
% (some pre-computation in the first section) 
% eigenvectors normalized and for stable representation
v = D*u;
v = v./vecnorm(v);

figure(3)
plot(diag(lambda),'x')
grid
title('eigenvalues (of the normalized Laplacian)')

F = 3;
figure(3)
plot(v(:,F),'x')

% simple approch: create two clusters with the Fiedler's vector
figure(4)
C = sign(v(:,F));
scatter(v(C == 1,F), v(C == 1,F+1), 'b', '*')
hold on;
scatter(v(C == -1,F), v(C == -1,F+1), 'r', '*')
hold off
grid
title('Communities')

% order nodes according to the value in the Fiedler's vector
[v1s,pos] = sort(v(:,F), 'descend');
Au1 = A(pos,pos);

% evaluate the conductance measure
a = sum(triu(Au1));
b = sum(tril(Au1)); %c = sum(triu(Au1),2)';
d = a+b;
D = sum(d);
assoc = cumsum(d); %assoc(S)
assoc = min(assoc,D-assoc); %min(assoc(S),assoc(S^c))
cut = cumsum(b-a);
conduct = cut./assoc;
conduct = conduct(1:end-1);

figure(5)
plot(conduct,'x-')
grid
title('conductance')
ylabel('conductance')

% identify the minimum -> threshold
[~,mpos] = min(conduct);
threshold = mean(v1s(mpos:mpos+1));
disp(['Minimum conductance: ' num2str(conduct(mpos))])
disp(['   Cheeger''s upper bound: ' num2str(sqrt(2*lambda(2,2)))])
disp(['   # of links: ' num2str(D/2)])
disp(['   Cut value: ' num2str(cut(mpos))])
disp(['   Assoc value: ' num2str(assoc(mpos))])
disp(['   Community size #1: ' num2str(mpos)])
disp(['   Community size #2: ' num2str(N-mpos)])

% show the result
figure(7)
plot(v(:,F),v(:,F+1),'.')
hold on;
plot(threshold*[1,1],ylim,'r-')
[~,mpos2] = min(conduct(1:100));
threshold2 = mean(v1s(mpos2:mpos2+1));
plot(threshold2*[1,1],ylim,'b-')
[~,mpos3] = min(conduct(300:400));
threshold3 = mean(v1s(mpos3:mpos3+1));
plot(threshold3*[1,1],ylim,'g-')
hold off
grid
title('communities')
legend('Nodes','Absolut minima','First minima','Center minima')




%% %%%%%%%%%%%%%%%%%%% PageRank-nibble approach %%%%%%%%%%%%%%%%%%%%%%%%%%%

% few useful things
d = sum(A); % degree vector
D = sum(d); % degrees sum
I = spdiags(ones(N,1),0,N,N); % identity matrix
Di = spdiags(1./sqrt(d'),0,N,N); % diagonal degrees square-rooted
L = I - Di*A*Di; % normalized Laplacian
M = A*Di*Di; % normalized adjacancy matrix

if mpos<N-mpos  % select seed node from the smaller group
    i = pos(1); % we select the more relevant from the perspective of the spectral approach
else
    i = pos(end);
end

% different start
%i = find(strcmp(names_of_MEPs,'BAY'));


q = zeros(N,1);
q(i) = 1; % teleport vector
c = 0.85;
r = (I-c*M)\((1-c)*q); % ranking vector
ep = 1e-3; % precision

% run PageRank-nibble
u = zeros(N,1); % starting point
v = q; % starting point
th = full(ep*d/D)'; % thresholds
count = 0; % exit counter
complexity = 0; % complexity value (# of operations)
ii = i; % starting index used for Push operation
while (count<N)
    if v(ii)>th(ii) % push if above threshold
        tmp = v(ii);
        u(ii) = u(ii)+(1-c)*tmp;
        v(ii) = 0;
        v = v + c*M(:,ii)*tmp;    
        complexity = complexity + d(ii); % update complexity
        count = 0; % reset the exit counter
    else % go to next entry if below threshold
        count = count + 1; % increase the exit counter
        ii = mod(ii,N)+1; % update the index used for Push
    end
end

% sweep wrt the ordering identified by v1
% reorder the adjacency matrix
[u1s,pos2] = sort(u,'descend');
Nmax = find(u1s>0,1,'last'); % discard nodes with 0 values (never used in Push)
Au1 = A(pos2,pos2(1:Nmax));
% evaluate the conductance measure
a = sum(triu(Au1));
b = sum(tril(Au1));
assoc = cumsum(a+b);
assoc = min(assoc,D-assoc);
cut = cumsum(b-a);
conduct = cut./assoc;
conduct = conduct(1:Nmax-1); 
% identify the minimum -> threshold
[~,mpos2] = min(conduct);
threshold2 = mean(u1s(mpos2:mpos2+1));
disp('PageRank-nibble approach')
disp(['   complexity/D: ' num2str((complexity/D))])
disp(['   epsilon: ' num2str(ep)])
disp(['   prec: ' num2str(norm(r-u,1))])
disp(['   Minimum conductance: ' num2str(conduct(mpos2))])
disp(['   # of links: ' num2str(D/2)])
disp(['   Cut value: ' num2str(cut(mpos2))])
disp(['   Assoc value: ' num2str(assoc(mpos2))])
disp(['   Community size #1: ' num2str(mpos2)])
disp(['   Community size #2: ' num2str(N-mpos2)])

% show sweep choice
figure(3)
plot(conduct, '-rx')
grid
ylabel('conductance')
title('Sweep choice')

% show network with partition
figure(4)
H = plot(u,v1s,'k.');
hold on
plot(u(pos2(1:mpos2)),v1s(pos2(1:mpos2)),'go')
plot(threshold2*[1,1],ylim,'g-')
plot(xlim,threshold*[1,1],'r-')
hold off
grid
ylabel('Fiedler''s eigenvector value')
xlabel('PageRank value')
title('Communities')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% K-MEANS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = diag(d.^-0.5);
L1 = spdiags(ones(N,1),0,N,N) - (D*A*D);

[u,lambda] = eigs(L1, N, 'smallestabs');

v = D*u;
v = v./vecnorm(v);

U = v(:,1:3);
% perform kmeans clustering on the matrix U
[IDX,~] = kmeans(U,3); 
% plot the eigen vector corresponding to the largest eigen value
%figure,plot(IDX)
figure
hold on;
for i=1:size(IDX,1)
    if IDX(i,1) == 1
        plot(v(i,2),v(i,3),'m+');
    elseif IDX(i,1) == 2
        plot(v(i,2),v(i,3),'g*');
    %elseif IDX(i,1) == 3
    %    plot(v(i,2),v(i,3),'b*');
    %elseif IDX(i,1) == 4
    %    plot(v(i,2),v(i,3),'g+');
    else
        plot(v(i,2),v(i,3),'b+');        
    end
end
hold off;
title('Clustering Results using K-means');
grid on;



%% plot with color from the clustering 
H = plot(GC, 'NodeColor','k','EdgeAlpha',0.1);
title('K-Means plot with color from the clustering');
%H.NodeLabel = index(pos);
%H.NodeLabel = names_of_MEPs;
p1 = find(IDX==1);
p2 = find(IDX==2);
p3 = find(IDX==3);
%p4 = find(IDX==4);
highlight(H,p1, 'NodeColor','g')
highlight(H,p2, 'NodeColor','r')
highlight(H,p3, 'NodeColor','k')
%highlight(H,p4, 'NodeColor','y')

disp(['Dim of cluster 1 = ' num2str(size(p1,1))])
disp(['Dim of cluster 2 = ' num2str(size(p2,1))])
disp(['Dim of cluster 3 = ' num2str(size(p3,1))])


%% plot with color from the political group
% three group: pro euro (1, 3, 4, 7, 8)
%              against euro (2, 5 )
%              no group (9) + (6) and (7) [that not have a clear position on
%              this topic]
% remember to uncomment at the beginning the line about group
figure(5)
H = plot(GC, 'NodeColor','k','EdgeAlpha',0.1);
title('K-Means plot with color from the groups');
%H.NodeLabel = names_of_MEPs;

% warning! The order of the calls of highlight change the results 
% because of the sovrapposition of some nodes 
for i=1:9
    p = find(group==i);
    if(i==9 || i==6  || i==7)
        highlight(H,p, 'NodeColor','c')
    elseif (i==2 || i==5)
        highlight(H,p, 'NodeColor','r')
    %elseif(i==1 || i==3 || i==4 || i==6 || i==7 || i==8)
    %    highlight(H,p, 'NodeColor','g')
    elseif(i==1 || i==3 || i==4 || i==8)
        highlight(H,p, 'NodeColor','g')
    end
end


%% support code to print out the clusters for other studies (SNA)
% fileID = fopen('p1.txt', 'w');
% 
% fprintf(fileID,'p1\n\n');
% fprintf(fileID,'%6.0f\n', p1);
% 
% fclose(fileID);

