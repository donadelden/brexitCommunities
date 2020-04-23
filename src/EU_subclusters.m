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

% plot graph
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

d = sum(A,2);
N = size(A,1);
D = diag(d.^-0.5);
L1 = spdiags(ones(N,1),0,N,N) - (D*A*D);

[u,lambda] = eigs(L1, N, 'smallestabs');

v = D*u;
v = v./vecnorm(v);

%% K MEANS
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

disp(['  #1: ' num2str(sum(IDX(:,1)==1))...
      '  #2: ' num2str(sum(IDX(:,1)==2))...
      '  #3: ' num2str(sum(IDX(:,1)==3))])
      

%% K MEANS on the biggest cluster
% please, before this run the section of kmeans on all the data
% this code is for three clusters
close all

% find biggest cluster
num_of_cl = 3;
new_cl_number = 4; % 5 max for the representation
p = zeros(num_of_cl,1);
for i=1:num_of_cl
    p(i) = sum(IDX==i);%size(find(IDX==i),1);
end

bc = find(p==max(p)); % index of biggest cluster
pos_bc = find(IDX==bc); % position of nodes of bc
A_bc = A(pos_bc,pos_bc); % generate the new adj matrix 

d_bc = sum(A_bc,2);
%A = W(pos,pos);
N_bc = size(A_bc,1);
D_bc = diag(d_bc.^-0.5);
L1_bc = spdiags(ones(N_bc,1),0,N_bc,N_bc) - (D_bc*A_bc*D_bc);
[u_bc,lambda_bc] = eigs(L1_bc, N_bc, 'smallestabs');

v_bc = D_bc*u_bc;
v_bc = v_bc./vecnorm(v_bc);

U_bc = v_bc(:,1:new_cl_number);
% perform kmeans clustering on the matrix U
[IDX_bc,~] = kmeans(U_bc,new_cl_number); 
figure(1)
hold on;
for i=1:size(IDX_bc,1)
    if IDX_bc(i,1) == 1
        plot(v_bc(i,2),v_bc(i,3),'m+');
    elseif IDX_bc(i,1) == 2
        plot(v_bc(i,2),v_bc(i,3),'g*');
    elseif IDX_bc(i,1) == 3
        plot(v_bc(i,2),v_bc(i,3),'b*');
    elseif IDX_bc(i,1) == 4
        plot(v_bc(i,2),v_bc(i,3),'c+');
    else
        plot(v_bc(i,2),v_bc(i,3),'m+');        
    end
end
hold off;
title('Clustering Results using K-means ON THE BIGGEST CLUSTER');
grid on;

%% visualize in graph
GC_bc = graph(A_bc);
figure(2)
H_bc = plot(GC_bc, 'NodeColor','k','EdgeAlpha',0.1);
%H_bc.NodeLabel = index(pos);
color = ['g', 'r', 'k', 'y', 'c'];
for i=1:new_cl_number
    p_bc = find(IDX_bc==i);
    highlight(H_bc, p_bc, 'NodeColor',color(i))
    disp(['Dim of cluster ' num2str(i) ' = ' num2str(size(p_bc,1))])
%     fileID = fopen(['p_bc_prebrexit' num2str(i) '.txt'], 'w');
%     
%     fprintf(fileID,['p_bc_' num2str(i) '\n\n']);
%     fprintf(fileID,'%6.0f\n', p_bc);
%     
%     fclose(fileID);
end