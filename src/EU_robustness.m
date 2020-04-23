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
H.NodeLabel = table2cell(names(pos,:));

L = numedges(GC);
nnode = numnodes(GC);
disp(['Node: ' num2str(nnode) ' - Edge: ' num2str(L)])
d = sum(A,1);



%% %%%%%%%%%%%%%%%%%%%% Random Removing of LINKS %%%%%%%%%%%%%%%%%%%%%%%%%%
% Warning! Due to the high number of links it's really slow, be patient!
sizesGC = [];
sizes = [];
N = size(A,1);
L = sum(sum(d))/2;
% random removal
Arem = A; % backup A matrix just in case...
% erase all the links is really slow, the last ones it's hard to find in a
% random way
number_of_link_to_remove = L; % L is the number of links

for count=1:number_of_link_to_remove
    % this way is good and fast until we've a big number of links
    if count<number_of_link_to_remove*0.7
        i = (floor(rand(1)* (N-1)) + 1);
        j = (floor(rand(1)* (N-1)) + 1);
        while(Arem(i,j) == 0)
            i = (floor(rand(1)* (N-1)) + 1);
            j = (floor(rand(1)* (N-1)) + 1);
        end
        Arem(i,j) = 0;
        Arem(j,i) = 0;
    else
    % this way is fast when we've a small amount of nodes because trying to
    % find a link in the whole matrix it's difficult, it's better to use
    % find function
        % find where are the ones in tre matrix
        ones_pos = find(Arem==1);
        % pick a random position
        i = (floor(rand(1)* (size(ones_pos,1)-1)) + 1);
        % delete the node
        Arem(ones_pos(i)) = 0;
        tmp = Arem';
        tmp(ones_pos(i)) = 0;
        Arem = tmp';
    end
    
    % in order to get better speed we check the size every 10 deletions
    if (mod(count,10)==0)
        G_tmp = graph(Arem);
        [~, bs] = conncomp(G_tmp);
        GC_size = max(bs);
        sizesGC(end+1) = GC_size;
        sizes(end+1) = count;
        disp(['rem: ' num2str(L - count) '  GC: ' num2str(sizesGC(end))...
            '   count: ' num2str(count)])
        
    end
end


figure(1)
plot(sizes/L, sizesGC)
xlabel('f (fraction of removed links)')
ylabel('size of the Giant Component')
title('Random removing of all the Links')
grid on

    
%% Random Removing of NODES
sizesGC2 = [];
sizes2 = [];
N = size(A,1);
% random removal
Arem = A>0; % backup A matrix and erase weigth just in case ...
number_of_node_to_remove = N-2;

for count=1:number_of_node_to_remove
    % select a random node
    i = (floor(rand(1)*(N-count-1))+1);
    % delete the node
    Arem(i,:) = []; 
    Arem(:,i) = [];
    % find the giant component
    G_tmp = graph(Arem);
    [~, bs] = conncomp(G_tmp);
    GC_size = max(bs);
    sizesGC2(end+1) = GC_size;
    sizes2(end+1) = N-size(Arem,1);
        
    
end

figure(1)
plot(sizes2/N, sizesGC2)
xlabel('f (fraction of removed nodes)')
ylabel('size of the Giant Component')
title('Random removing of Nodes')
grid on

%% Attack to the Nodes
% remove nodes starting from the hubs (decreasing degree)
sizesGC3 = [];
sizes3 = [];
N = size(A,1);
Arem = A>0; % backup A matrix and erase weigth just in case ...
number_of_node_to_remove = N-2;

for count=1:number_of_node_to_remove
    d = sum(Arem);
    % select node with higher degree
    i = find(d==max(d));
    i = i(1); % just in case there is more of one with higher degree
    % delete the node
    Arem(i,:) = []; 
    Arem(:,i) = [];
    
    % find the giant component    
    G_tmp = graph(Arem);
    [~, bs] = conncomp(G_tmp);
    GC_size = max(bs);
    sizesGC3(end+1) = GC_size;
    sizes3(end+1) = N-size(Arem,1);
end

figure(1)
plot(sizes3/N, sizesGC3, '-g')
%hold on
%plot(sizes2/N, sizesGC2, 'r-')
xlabel('f')
ylabel('size of the Giant Component')
title('Attack/random removes on nodes')
%legend('Attack', 'Random removes');
grid on
