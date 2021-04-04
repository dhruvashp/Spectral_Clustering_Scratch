% A single spectral clustering run

X = readmatrix('X_2.csv');             
y_true = readmatrix('s_2.csv');

% test is done based on aca2 data, as such code and everything else stays
% the same for both aca2 and aca5 so it doesn't matter which is used


gamma = 0.5;             % fixed values for test
k = 5;                   % fixed values for test


% For Kernel K

% Note : X has data points along columns, for X an m x n matrix, we have n
% datapoints with dimension of each being m

% Kernel K will thus be n x n

K = zeros(size(X,2),size(X,2));           % initializing

for i = 1:size(X,2)
    for j = 1:size(X,2)
        K(i,j) = exp((-gamma)*((norm(X(:,i)-X(:,j)))^2));
    end
end

% The Kernel has been built
% Now to build the Weight Matrix

W = K;              % we will also save K, in case we need it later

for j = 1:size(W,2)
    
    [S,I] = sort(W(:,j),'descend');
    
    for i = 1:k
        
        W(I(i),j) = S(i);           % this step is not really needed, but done for validation (since W already has the top k values in its original index)
        
    end
    
    for i = k+1:size(I,1)
        
        W(I(i),j) = 0;              % this is the only step actually needed, as it makes those indices of W (which has K) 0 that don't contain the top k values in that column
        
    end
    
    
    
end


% the output has been tested and verified using smaller matrices in Test


A = (W + transpose(W))/2;         % symmetrizing W



% k : Columnar top k entries of K to be selected. Note k = 1 will yield a
% diagonal W, since the diagonal entries are always 1 in K (largest). Thus
% we interpret k given in HW to be k including this diagonal 1, which, will
% always be redundantly included. Thus this k here is inclusive of this
% 'ever-present' diagonal 1. We can define k' which is all non-redundant
% top entries in each column. k' = 1, then, is defined, and is equivalent
% to k = 2. Thus k = k' + 1, k starts from 2 and k' from 1. Here we've
% taken k to be the redundant k, since in HW, the k mentioned does start
% from 2.

% K : original kernel, with 1's along the diagonal. In the paper it is 0,
% but we assume 0/1 diagonal entries both valid. Here, thus, the diagonal
% entries are 1.

% W : sparse kernel of K, based on k, non-symmetric, always has 1's on the
% diagonal

% A : final affinity, weight, similarity matrix, symmetric, sparse



% ASSUMPTIONS :
% We assume the diagonals being 1 and the 'k-top' selection, is in line
% with the details provided in the paper. 




% Making D

D = zeros(size(X,2),size(X,2));

for i = 1:size(X,2)
    D(i,i) = sum(A(i,:));
end

% Making L                  

L = (D^(-1/2))*A*(D^(-1/2)); 

% This is the L mentioned in the paper. Difference between :
% L   (paper)             : We select largest eigenvectors
% I-L (class)             : We select smallest eigenvectors

% Selecting smallest c (clusters) eigenvalues of I-L corresponds to
% selecting largest c (clusters) eigenvalues of L

% I-L was symmetric and PSD, L is symmetric. 
% If eigenvalues of I-L are [0,1], only then, L is PSD
% Thus if I-L has even a single eigenvalue larger than 1, than, L won't be
% PSD.

% When I-L eigenvalues are arranged in ascending order, we start from the
% smallest eigenvalue, and, keep appending corresponding eigenvectors till
% we see a sudden rise in eigenvalue

% Similarly, when L eigenvalues, are also arranged in ascending order, we
% start from the largest eigenvalue, append corresponding eigenvectors of
% these large eigenvalues till we see a sudden drop in eigenvalue



% I-L (class) : Start from bottom. Check rise. Decide cluster number. 
% L   (paper) : Start from top. Check drop. Decide cluster number.
%       (after arranging all eigenvalues in ascending order)


[Z,lam] = eig(L);

% MATLAB does not generate sorted Z and lam for Diagonalization like it
% does for SVD. Thus we need to sort Z and lam in a corresponding fashion

[dum,indx] = sort(diag(lam),'descend');

lam_s = lam(indx,indx);

Z_s = Z(:,indx);

% lam_s and Z_s have been arranged, in a descending and corresponding
% fashion

% NOTE : As such, for a generic clustering problem whereby no labels are
% given, obviously, the cluster number c is to be decided based on the
% 'maximal dip or rise' analysis. Here, since, we know the labels, if we
% assume label distinctness to correspond to 'clustered entities' that are
% dissimilar, we may obviously obtain the best results from clustering if
% we cluster with the cluster number = maximum number of distinct labels

% Thus, 

% clusters = maximum number of distinct labels in our original label set
% (y_true)



% Selecting Cluster Number 

clus = size(unique(y_true),2);


% unique returns a vector that contains all the classes only once, size
% calculates the length of such a vector, thus calculating the maximum
% number of distinct classes, obviously, without repeating the same class
% twice


eigvec_selec = Z_s(:,1:clus);
eigval_selec = lam_s(1:clus,1:clus);

data = eigvec_selec;          % copying selected eigenvectors to data

data = real(data);

% normalizing eigvec_selec along its rows, before clustering

for i = 1:size(data,1)

    if norm(data(i,:)) == 0
        continue
    else
        data(i,:) = data(i,:)/norm(data(i,:));
    end

end



% data contains the final data on which we have to perform k-means
% clustering (data points are in the rows)

% each row contains a data point vector to cluster



% K-Means on data   (with clusters = clus)      


% NOTE : The paper uses a K-Means algorithm that selects cluster centroids
% orthogonally from one another, which, given the discussions in class,
% makes sense.

% The difference majorly lies in the K-Means initiation. According to
% paper, initially centroids are assigned. The first centroid is picked
% randomly from the data points itself (the data rows). All the subsequent
% centroids, in the initiation, are also picked from the data points, such
% that all centroids are as orthogonal to each other as possible.
% Subsequent procedure stays same as in traditional K-Means.

% Here, however we'll use the traditional K-Means algorithm, as has been
% defined, in general. Since the traditional K-Means is obviously more
% widely applicable, applying it will be more utilitarian in general.








% TRADITIONAL K-MEANS ON DATA (data)


epsilon = 0.5;


ini_cluster_occur_count = zeros(clus,1);                      % column vector, row 1 = count of cluster 1 occurence,..........row clus = count of cluster clus occurence

while sum((ini_cluster_occur_count == 0)) ~= 0                % this while loop ensures that the initial random cluster assignment is such that each cluster is assigned 
                                                              % to at least a data point, and that there is no such cluster that isn't assigned to any data point

    
    ini_cluster_assign = randi([1 clus],1,size(data,1));      % initial clusters are assigned randomly, in traditional K-Means
                                                              % row vector

    [ini_cluster_assign_sort,ini_cluster_assign_indx] = sort(ini_cluster_assign,'ascend');


    term = 1;
    for m = 1:clus


        count = 0;

        for n = term:size(ini_cluster_assign_sort,2)

            if ini_cluster_assign_sort(1,n) == m
                count = count + 1;
            else
                term = n;
                break
            end

        end

        ini_cluster_occur_count(m,1) = count;

        if sum(ini_cluster_occur_count) == size(ini_cluster_assign_sort,2)
            break
        end

    end
    
    
    
    
    
    
end


cluster_assign = ini_cluster_assign;                           % this is the actual initial random cluster assignment, which ensures that each cluster appears at least once 
cluster_assign_sort = ini_cluster_assign_sort;                 % the above random initial assignment sorted
cluster_assign_indx = ini_cluster_assign_indx;                 % the original index for relevant data segregation
cluster_occur_count = ini_cluster_occur_count;                 % the occurence count of each cluster, ordered, column vector


% the variables in the final while loop trial are the same variables which
% are stored, since, obviously, those correct variables ensured the while
% loop termination


% NOTE : The while loop, may, seem trivial, especially when number of
% clusters are small compared to number of data points, since, it then
% would be rather rare for all clusters to not be assigned at least once

% However for generality, we've included this possibility in our algorithm


cluster_centroids = zeros(clus,clus);                          % cluster centroids, row 1 = centroid of cluster 1,...........row clus = centroid of cluster clus
                                                               % initializing
                                                               
pre = 1;

for t = 1:clus
    
    cluster_centroids(t,:) =  mean(data(cluster_assign_indx(1,pre:pre+cluster_occur_count(t,1)-1),:),1);
    pre = pre+cluster_occur_count(t,1);
    
end


% loss function calculation (optimizing quantity)

prece = 1;
loss = 0;

for q = 1:clus
    
    loss = sum(pdist2(data(cluster_assign_indx(1,prece:prece+cluster_occur_count(q,1)-1),:),cluster_centroids(q,:),'euclidean'))+loss;
    prece = prece+cluster_occur_count(q,1);
    
end


%disp(loss)                % checking


% the required loss function to monitor K-Means iterations has been
% obtained


% the above code, generates the cluster_centroids matrix, ordered

% to do so, we slice the indx vector at positions that give the
% corresponding clusters data indices

% using the sliced indices, we slice data matrix, and then obtain the mean
% over the row axis (axis = 1) to compute the cluster corresponding
% centroid




% classifying a data point in data, now, to a new cluster

% each data points distance will be calculated to all the cluster centers
% (centroids) and the data point will be assigned a new cluster, based on
% which clusters cluster center is closest to that data point



% looping over steps, down below, to keep on updating the clusters till
% loss function stabilizes



loss_pre = 0;




while abs(loss - loss_pre) > epsilon    
    
    loss_pre = loss;
    
    cluster_assign_updated = zeros(1,size(data,1));

    for l = 1:size(data,1)

        cluster_distance = zeros(clus,1);

        for p = 1:clus
            cluster_distance(p,1) = norm(data(l,:)-cluster_centroids(p,:));
        end

        [clus_sel_min_dist,clus_sel] = min(cluster_distance); 

        cluster_assign_updated(1,l) = clus_sel;

    end
    
    
    cluster_assign = cluster_assign_updated;
    
    
    
    
    
    
    
    
    cluster_occur_count = zeros(clus,1);
    
    [cluster_assign_sort,cluster_assign_indx] = sort(cluster_assign,'ascend');


    term = 1;
    for m = 1:clus


        count = 0;

        for n = term:size(cluster_assign_sort,2)

            if cluster_assign_sort(1,n) == m
                count = count + 1;
            else
                term = n;
                break
            end

        end

        cluster_occur_count(m,1) = count;

        if sum(cluster_occur_count) == size(cluster_assign_sort,2)
            break
        end

    end
    
    cluster_centroids = zeros(clus,clus);
    
    pre = 1;

    for t = 1:clus

        cluster_centroids(t,:) =  mean(data(cluster_assign_indx(1,pre:pre+cluster_occur_count(t,1)-1),:),1);
        pre = pre+cluster_occur_count(t,1);

    end

    prece = 1;
    loss = 0;

    for q = 1:clus

        loss = sum(pdist2(data(cluster_assign_indx(1,prece:prece+cluster_occur_count(q,1)-1),:),cluster_centroids(q,:),'euclidean'))+loss;
        prece = prece+cluster_occur_count(q,1);

    end
    
    
    % loss is updated, based on the cluster assignment which is updated
    % initially
    
    
    
    %disp(loss)         % checking
    
    
    

end


% cluster_assign contains the final clusters after loss stabilization

% loss contains the final loss after loss stabilization

% cluster_centroids contains the corresponding cluster centroids

% cluster_assign_sort and cluster_assign_indx contain the corresponding
% cluster_assign, sorted, and with original index

% cluster_occur_count contains cluster count for each cluster


% NOTE : Only a single run or 'restart' performed to save computational
% time


% Adding more than one restart will hurt both the memory and the
% computation






% Assigning classes to clusters

% Each cluster will now be polled to assign a class to a cluster

% The class (true class) of all the data points in each cluster will be
% polled and a class assigned to a cluster




cluster_class = zeros(clus,clus);          % row contains corresponding cluster, each column contains count of classes in each cluster

pri = 1;
for u = 1:clus
    slice_indx = cluster_assign_indx(1,pri:pri+cluster_occur_count(u,1)-1);
    y_true_slice = y_true(1,slice_indx);
    pri = pri+cluster_occur_count(u,1);
    for cls = 1:clus
        loc = 0;
        for len = 1:size(y_true_slice,2)
            if y_true_slice(1,len) == cls
                loc = loc + 1;
            end
        end
        cluster_class(u,cls) = loc;
    end
end


%disp(cluster_class)        % checking

[redundant,clus_cls_vec] = max(cluster_class,[],2);    % clus_cls_vec contains the relevant vector assigning class to each cluster

%disp(clus_cls_vec)          % checking


y_predicted = zeros(1,size(y_true,2));

for g = 1:size(cluster_assign,2)
    y_predicted(1,g) = clus_cls_vec(cluster_assign(1,g),1);
end

% y_predicted generated


mis = 0;
for w = 1:size(y_true,2)
    if y_predicted(1,w) ~= y_true(1,w)
        mis = mis + 1;
    end
end

err = (mis/size(y_true,2))*100;

disp('The percentage error for this test case is : ')
disp(err)





