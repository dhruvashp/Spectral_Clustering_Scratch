% Using inbuilt K-Means provided by MATLAB

% NOTE :
% All comments have been removed for simplicity
% Inbuilt K-Means instead of scratch used here (with 3 restarts)





X = readmatrix('X_2.csv');             
y_true = readmatrix('s_2.csv');




gamma = 0.5;             
k = 5;                   




K = zeros(size(X,2),size(X,2));          

for i = 1:size(X,2)
    for j = 1:size(X,2)
        K(i,j) = exp((-gamma)*((norm(X(:,i)-X(:,j)))^2));
    end
end


W = K;              

for j = 1:size(W,2)
    
    [S,I] = sort(W(:,j),'descend');
    
    for i = 1:k
        
        W(I(i),j) = S(i);           
        
    end
    
    for i = k+1:size(I,1)
        
        W(I(i),j) = 0;              
        
    end
    
    
    
end





A = (W + transpose(W))/2;         













D = zeros(size(X,2),size(X,2));

for i = 1:size(X,2)
    D(i,i) = sum(A(i,:));
end

              

L = (D^(-1/2))*A*(D^(-1/2)); 







[Z,lam] = eig(L);


[dum,indx] = sort(diag(lam),'descend');

lam_s = lam(indx,indx);

Z_s = Z(:,indx);






clus = size(unique(y_true),2);





eigvec_selec = Z_s(:,1:clus);
eigval_selec = lam_s(1:clus,1:clus);

data = eigvec_selec;         

data = real(data);

for i = 1:size(data,1)
    
    if norm(data(i,:)) == 0
        continue
    else
        data(i,:) = data(i,:)/norm(data(i,:));
    end
    
    

end








% K-Means       (using inbuilt MATLAB function)





cluster_assign_int = kmeans(data,clus,'Replicates',3,'Start','uniform');               % performed thrice (restarted thrice), each time the 
                                                                                       % initial cluster centroids picked at random from data
                                                                                       
cluster_assign = transpose(cluster_assign_int);

[cluster_assign_sort,cluster_assign_indx] = sort(cluster_assign,'ascend');


cluster_occur_count = zeros(clus,1);

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





cluster_class = zeros(clus,clus);
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




[redundant,clus_cls_vec] = max(cluster_class,[],2);    




y_predicted = zeros(1,size(y_true,2));

for g = 1:size(cluster_assign,2)
    y_predicted(1,g) = clus_cls_vec(cluster_assign(1,g),1);
end



mis = 0;
for w = 1:size(y_true,2)
    if y_predicted(1,w) ~= y_true(1,w)
        mis = mis + 1;
    end
end

err = (mis/size(y_true,2))*100;

disp('The percentage error for this test case is : ')
disp(err)





