% Using an SVM model

clear all
close all
clc
 
%rng default

% 1) Prepare dataset

% load('/Volumes/WD Ezra/Data/Synaptosomes/Experiment_37C/Results/results_combined_after_overlap_threshold.mat');
% load('/Volumes/WD Ezra/Data/Synaptosomes/Experiment_37C/Results/Results_combined/ripley/clustersizes.mat');
% load('/Volumes/WD Ezra/Data/Synaptosomes/Experiment_37C/Results/Results_combined/ripley/interclusterdist.mat');
% load('/Volumes/WD Ezra/Data/Synaptosomes/Experiment_4C/Results/results_combined_after_overlap_threshold.mat');
% load('/Volumes/WD Ezra/Data/Synaptosomes/Experiment_4C/Results/Results_combined/ripley/clustersizes.mat');
% load('/Volumes/WD Ezra/Data/Synaptosomes/Experiment_4C/Results/Results_combined/ripley/interclusterdist.mat');
% load('/Volumes/WD Ezra/Dump/results_combined_after_overlap_threshold.mat');
load('/Users/Ezra/Desktop/data_37C_classified.mat');


% Remove rows that have -1 value in any of the Ripley's columns
data(ismember(data.clustersizeRC,-1),:)=[];
data(ismember(data.clustersizeGC,-1),:)=[];
data(ismember(data.clustersizeBC,-1),:)=[];
data(ismember(data.interclusterdistRG,-1),:)=[];
data(ismember(data.interclusterdistRB,-1),:)=[];
data(ismember(data.interclusterdistGB,-1),:)=[];


y = cell2mat(table2array(data(:,3)));

%X = [data(:,1:7) data(:,12:13) data(:,8:11) data(:,14:end)];

%X = [array2table(grp2idx(table2array(data(:,1))),'VariableNames',{'condition'}) data(:,7:end)]; % only select usefull measurements from data table
%X = [data(:,7:end)]; % only select usefull measurements from data table
X = [array2table(grp2idx(table2array(data(:,1))),'VariableNames',{'condition'}) data(:,12:end)]; % only select usefull measurements from data table

%%
% 2) Divide in training and test set (e.g. 80/20 for train/test ratio)

% Randomly permutate integers between 1 and length(y)
rand_num = randperm(length(y));

% Randomly divide into train and test data
X_train = X(rand_num(1:round(0.80*size(data,1))),:);
y_train = y(rand_num(1:round(0.80*size(data,1))),:);
X_test  = X(rand_num(round(0.80*size(data,1))+1:end),:);
y_test  = y(rand_num(round(0.80*size(data,1))+1:end),:);

% % Divide into train and test data
% X_train = X(1:80,:);
% y_train = y(1:80,:);
% X_test  = X(81:end,:);
% y_test  = y(81:end,:);

%%
% 3) Prepare validation set out of training set (k-fold cross-validation)
c = cvpartition(y_train,'k',5);

%%
% 4) Feature selection

% Some display options
opts = statset('display','iter');

% Fun compares the prediction of the labels of some test data feeded into 
% the model generated by fitcsvm with the true test labels. So it will
% check how 'wrong' the model was by summing the errors
fun = @(train_data, train_labels, test_data, test_labels)...
    sum(predict(fitcsvm(train_data, train_labels, 'KernelFunction','rbf'), test_data) ~= test_labels);

% Perform the feature selection (This one uses forward sequential feature
% selection by default. You can add 'direction','backward' as options;)
[fs,history] = sequentialfs(fun, table2array(X_train), y_train, 'cv', c, 'options', opts, 'nfeatures', 3);
% You can also specify which columns are definitely included, or definitely
% not included.
% 'fs' gives you which final columns are included
% 'history' gives you the history of different combinations of columns that
% were tried out before

%%
% 5) Finding the best hyper parameters for the classification

% Use only the features that were found to be most useful for
% classification
X_train_with_best_features = X_train(:,fs);

% Get the SVM model
SVMModel = fitcsvm(X_train_with_best_features,y_train);


% 6) Test the data with the test set
X_test_with_best_features = X_test(:,fs);
accuracy = sum(predict(SVMModel, X_test_with_best_features) == y_test)/length(y_test)*100;
disp(['Model predicted the right label ' num2str(accuracy) '% of the time.'])

Y_test = predict(SVMModel, X_test_with_best_features);

%% 7) Show the hyperplane in the model

x1 = table2array(X_train_with_best_features(y_train==0,:));
x2 = table2array(X_train_with_best_features(y_train==1,:));

figure(1);
scatter3(x1(:,1),x1(:,2),x1(:,3),'r');
hold on
scatter3(x2(:,1),x2(:,2),x2(:,3),'b');
hold off
xlabel(SVMModel.PredictorNames(1))
ylabel(SVMModel.PredictorNames(2))
zlabel(SVMModel.PredictorNames(3))


x3 = table2array(X_test_with_best_features(Y_test==0,:));
x4 = table2array(X_test_with_best_features(Y_test==1,:));

figure(2)
scatter3(x3(:,1),x3(:,2),x3(:,3),'r');
hold on
scatter3(x4(:,1),x4(:,2),x4(:,3),'b');
hold off
xlabel(SVMModel.PredictorNames(1))
ylabel(SVMModel.PredictorNames(2))
zlabel(SVMModel.PredictorNames(3))


% figure;
% hgscatter = gscatter(table2array(X_train_with_best_features(:,1)),...
%                      table2array(X_train_with_best_features(:,2)),...
%                      y_train);
% hold on;
% h_sv = plot(SVMModel.SupportVectors(:,1),...
%             SVMModel.SupportVectors(:,2),...
%             'ko','markersize',8);
% xlabel(SVMModel.PredictorNames(1))
% ylabel(SVMModel.PredictorNames(2))
% set(gca,'fontsize',14);
% legend('compact','disperse')