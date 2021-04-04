% Data Preparation
% As has been described in HW, we will prepare the data, which is to be
% clustered
% Code will be written for X and s (data and true class) in a generic
% fashion

% Then this code will be RUN once, after importing aca2, and again after
% importing aca5 (X and s in both)

% This final data to be clustered will then be stored before being utilized
% again for subsequent portions of the HW


% Normalizing Columns of X (as data is in columns)

for i = 1:size(X,2)
    X(:,i) = X(:,i)/norm(X(:,i)); %#ok<SAGROW>
end

% X is now normalized

X = X(:,1:2:end);
s = s(1:2:end);

writematrix(X,'X_aca5.csv')          % first X_aca2 was done, then X_aca5 was done
writematrix(s,'s_aca5.csv')          % first s_aca2 was done, then s_aca5 was done