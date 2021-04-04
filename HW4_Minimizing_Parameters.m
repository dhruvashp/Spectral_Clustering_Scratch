e_2 = readmatrix('Error_Matrix_aca2');
e_5 = readmatrix('Error_Matrix_aca5');

gamma = [[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9],1:100];
k = 2:50;

[m_e_2,i_e_2] = min(e_2,[],1);

[ovr_min_err_2,col_e_2] = min(m_e_2,[],2);

row_e_2 = i_e_2(1,col_e_2);

gamma_2 = gamma(1,row_e_2);
k_2 = k(1,col_e_2);


[m_e_5,i_e_5] = min(e_5,[],1);

[ovr_min_err_5,col_e_5] = min(m_e_5,[],2);

row_e_5 = i_e_5(1,col_e_5);

gamma_5 = gamma(1,row_e_5);
k_5 = k(1,col_e_5);



disp('Results for aca2 : ')
disp('Optimal gamma, aca2 : ')
disp(gamma_2)
disp('Optimal k, aca2 : ')
disp(k_2)
disp('Corresponding minimum % error, aca2 : ')
disp(ovr_min_err_2)




disp('Results for aca5 : ')
disp('Optimal gamma, aca5 : ')
disp(gamma_5)
disp('Optimal k, aca5 : ')
disp(k_5)
disp('Corresponding minimum % error, aca5 : ')
disp(ovr_min_err_5)









