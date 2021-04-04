There are a total of three test files that were executed prior to
executing the final loop over gamma and k

All three were run over aca2 data


gamma = 0.5 (for all three)
k = 5 (for all three)




Scratch_1 and Scratch_2



These files implement, over a single gamma and k loop (with aforementioned values), the Spectral Clustering algorithm, with the K-Means implemented from scratch without using the inbuilt function in MATLAB. 

The K-Means here, implemented from scratch, was restarted only once, however, additional code to restart it multiple times (for global convergence) would take only 2-3 lines of additional code to loop over the single restart K-Means, storing clusters in each restart with the loss in each restart.

Difference between Scratch_1 and Scratch_2

Scratch_1 : The cluster initiation is done by assigning clusters randomly to each data point, ensuring, that all clusters are assigned at least once

Scratch_2 : The cluster initiation is done by selecting data points randomly as cluster centroids (obviously, no data point is repeated in the cluster centroids)








Inbuilt

As the name suggests, this uses the inbuilt MATLAB function for K-Means. The test case has been restarted thrice, internally.









NOTE : 
For gamma = 0.5 and k = 5, the code written from scratch (for K-Means) performed just as well as the code where inbuilt MATLAB function was used. Results could have further been improved for the scratch code if it was also restarted. Again that modification would have been rather minor, but, since the K-Means algorithm gist was captured in only a 'single restart', I chose not to perform multiple restarts for the code from scratch.
