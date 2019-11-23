% % % %  Supervised Classification of 5435 Audio Signals    % % % %
% % % %      Christian Grenier, Tianyi Zhou, Zach Neveu     % % % %
% % % %                    04/15/2019                       % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% This code reads the 5435 audio files and their labels, runs PCA
%   on extracted 107 features, performs 5 supervised classification 
%   with 10-fold cross validation and 80-20 holdout validation. 
%   Classification results are evaluated using confusion matrix.
%
%   
%   Tip: To save time, you could skip Pre-processing and Feature 
%           Extraction section, and load normFeatures.mat at the
%           beginning in PCA section. Then run the classification
%           section to see the results.
%% Clean workspace
clear
clc
close all

%% Pre-processing
fprintf('\n.......Loading and sorting data......\n')
tic % Start the timer
% ---------------------------------
% Import audio data from folder 
% ---------------------------------

% Codes written on WINDOWS. Path should be changed under other OS.
dataDir = 'RAWDATAunzipped\train\Train';
dataStore = audioDatastore(dataDir,'FileExtensions','.wav');


% --------------------------------------------------
% Imported data is in the order 0,1,10,100...
% We need to sort files in natural order 0,1,2,3...
% Sorted result saved in dataStore(audioDatastore)
% --------------------------------------------------

% Grab the number in filename and save it to ID(double array)
dataFiles = cell2table(dataStore.Files,'VariableNames',{'FilePath'});
NUMdata = length(dataFiles.FilePath);
ID = zeros(NUMdata,1);
for i = 1:NUMdata 
    [~,fn,~] = fileparts(char(dataFiles.FilePath(i)));
    ID(i) = str2num(fn);
end

% Sort files in dataStore according to ID in ascending order
dataFiles = [table(ID) dataFiles];
[dataFiles,changedIdx] = sortrows(dataFiles);
dataStore.Files = dataStore.Files(changedIdx); 

labelFile = readtable('RAWDATAunzipped\train\train.csv');
dataStore.Labels = categorical(labelFile.Class);

% ----------------------------------------------------------
% Randomly select audio signal for each class and visualize
% ----------------------------------------------------------

% Find 10 class names and get idx for first sample in that class
[classNames, idx] = unique(dataStore.Labels);

for i=1:10
    [signal,fs] = audioread(dataStore.Files{idx(i)});
    nChannel = numel(signal)/fs;    % Find number of channles
    dt = 1/fs;
    t = 0:dt:(length(signal)*dt)-dt;
    
    % Plot the audio signal
    figure(i)
    if nChannel==2
        plot(t,signal(:,1),t,signal(:,2))   % Plot both channels
    else
        plot(t,signal)                      % Plot single channel
    end
    
    class = strrep(string(classNames(i)), '_', ' ');    % Change underscore to space in class name
    title(class)
    xlabel('Seconds')
    ylabel('Amplitude')
    legend('Channel 1','Channel 2')
    axis([0 4 -1 1])                        % Set axis limits
    set(gcf,'position',[400,400,800,200])   % Set size of the figure
    saveas(figure(i),strcat('Signal_',class),'png')
end


toc % Stop the timer and display elapsed time
%% Feature Extraction
% ------------------------------------------------------------------------
% We use 107 features in total to represent one audio signal
%   2 for basic info:     Sample Rate, Duration;
%   2 for loudness:       Integrated Loudness, Loudness Range
%   2 for pitch:          Average Pitch, Average Change in Pitch
%   3 for Harmonic Ratio: Average HR, Variance in HR, Average Change in HR
%  98 for MFCC:           7 time windows for each audio signal, 1 log
%                            energy value and 13 mel frequency cepstral
%                            coefficients for each window
% See below for more details.
% ------------------------------------------------------------------------

% Create empty arrays to store results
SampleRate = zeros(NUMdata,1);
AudioFiles = cell(NUMdata,1);
Duration = zeros(NUMdata,1);

Loudness = zeros(NUMdata,1);
LdnRange = zeros(NUMdata,1);

Pitch_avg = zeros(NUMdata,1);
Pitch_diff_avg = zeros(NUMdata,1);

HR_avg = zeros(NUMdata,1);
HR_var = zeros(NUMdata,1);
HR_diff_avg = zeros(NUMdata,1);

padding = zeros(NUMdata,1);
AudioFilesPad = cell(NUMdata,1);
MFCC = cell(NUMdata,1);
MFCCmat = zeros(NUMdata,98); % N_Window * N_coeff = 7*(13+1) = 98

fprintf('\n......Reading audio files and extracing features......\n')
tic
for i=1:NUMdata
    % Read audio file from dataStore
    [audioIn, Info] = read(dataStore);
    audioIn = audioIn(:,1);                % Dump 2nd channel of audio data
    
    %   We ignore the second channel of audio data to reduce complexity,
    %   but feel that it would improve the accuracy of our results if
    %   properly integrated.
    
    % -------------------------
    % Basic Info
	% -------------------------

    SR = Info.SampleRate;
    SampleRate(i) = SR;
    Duration(i) = length(audioIn)/SR;
    AudioFiles(i) = {audioIn};
    fprintf('Reading file %d. Completed %.1f%%.\n',...
        i,progress(dataStore)*100)

    % Loudness    
    % -------------------------------------------------------
    % integratedLoudness is a function in Audio Toolbox, it
    %   returns integrated loudness as a scalar value and 
    %   loudness range, a measurement of change in loudness 
    %   across the audio file, max - min, as a scalar value
    % -------------------------------------------------------

    if Duration(i) < 0.4        % Ldn is returned empty if audio
        Loudness(i) = NaN;      %    is shorter than 0.4s
        LdnRange(i) = NaN;       
    elseif Duration(i) < 3
        Loudness(i) = integratedLoudness(audioIn,SR); 
        LdnRange(i) = NaN;      % LR is returned empty if audio 
    else                        %   is shorter than 3s
        [Loudness(i),LdnRange(i)] = integratedLoudness(audioIn,SR); 
    end

    % Pitch
    % -------------------------------------------------------
    % pitch is a function in Audio Toolbox, it returns the 
    %   fundamental frequency of audio signal over time. 
    %   Average Pitch measured the mean pitch across the audio.
    %   Average Change in Pitch measures how much pitch changes 
    %   over time.
    % -------------------------------------------------------
    if Duration(i) < 0.052          % Pitch is returned empty if audio
        Pitch = NaN;                %    is shorter than 0.052s
        Pitch_avg(i) = NaN;
        Pitch_diff_avg(i) = NaN;
    else 
        Pitch = pitch(audioIn,SR);
        Pitch_avg(i) = mean(Pitch);
        Pitch_diff_avg(i) = mean(abs(diff(Pitch)));
    end
    
    % Harmonic Ratio
    % -------------------------------------------------------
    % harmonicRatio (HR) is a function in Audio Toolbox, it returns 
    %   the ratio of the energy at the fundamental frequency of the 
    %   signal to the total energy of the audio. This is a windowed 
    %   measurement and the lenght of the feature depends on the 
    %   window length and overlap, but we distilled it down to 3 usable
    %   features, similar to pitch. The mean of the HR, the variance of
    %   the HR, and average magnitude of the rate of change to get an idea
    %   of how it changes over the signal.
    % -------------------------------------------------------
    HR = harmonicRatio(audioIn,SR);
    HR_avg(i) = mean(HR);
    HR_var(i) = var(HR);
    HR_diff_avg(i) = mean(abs(diff(HR)));
    
    % MFCC
    % -------------------------------------------------------
    % mfcc is a function in Audio Toolbox, it returns 
    %   mel frequency cepstral coefficients for an audio file
    %   over time. We padded the short audio to 4s, and added small SNR
    %   gaussian white noise so that the MFCC function doesn't return
    %   -inf when taking the log of the DFT. Then we manually set the 
    %   window length and the overlap length so that the number
    %   of windows is 7. For each time window, it calculates 
    %   1 log energy value and 13 mfcc.
    % -------------------------------------------------------
    
    % Padding short audio to prepare for MFCC extraction
    padding(i) = uint32((4 - Duration(i))*SR);
    padArray = zeros(padding(i),1);
    AudioFilesPad(i) = {awgn([audioIn; padArray],10e2)};
    
    WL = round(SR*1);   % Window Length ~1s
    OL = round(WL/2); % Overlap Length ~0.5s
    MFCC(i) = {mfcc(AudioFilesPad{i},SR,'WindowLength',WL,'OverlapLength',OL)};
    
    %figure(i)
    %pcolor(MFCC{i}')
    
    % Reshape 7x14 MFCC{i} array to 1x98 MFCCmat(i,:) array
    MFCCmat(i,:) = reshape(MFCC{i},[1,98]);

end
toc

% Construct a table containing all information
dataFiles = [dataFiles table(AudioFiles) table(SampleRate) ...
    table(Duration) table(Loudness) table(LdnRange) table(Pitch_avg) ...
    table(Pitch_diff_avg) table(HR_avg) table(HR_var) ...
    table(HR_diff_avg) table(MFCCmat)];

% Construct an array containing only features
Features = [SampleRate Duration Loudness LdnRange Pitch_avg Pitch_diff_avg ...
    HR_avg HR_var HR_diff_avg MFCCmat];

% Customized zscore function to deal with missing data in feature array
zscore_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')),...
    std(x, 'omitnan'));
normFeatures = zscore_xnan(Features);

save normFeatures.mat normFeatures

%% Dimensionality Reduction - PCA

% Uncomment this line if you want to use pre-calculated normFeatures
load normFeatures.mat          

rng(2);         % For reproducibility

fprintf('\n......Starting PCA......\n')
tic

% Use ALS algorithm to deal with missing data
[~,score,latent] = pca(normFeatures,'Algorithm','ALS'); 

% Find number of principal components
tot = sum(latent);
perc = cumsum(latent)./tot;
NumOfPC = find(perc>0.95, 1 );  % Take components contributing over 95%

% Plot the Pareto Chart for all components
figure(11)
bar(latent)
hold on
yyaxis right
plot(1:length(latent),perc,'o-')
hold off
title('Pareto Chart for Principal Components')
saveas(figure(11),'Pareto Chart','png')

% Find final feature table after dimensionality reduction
Data_array = score(:,1:NumOfPC);        % data without labels

labelFile = readtable('RAWDATAunzipped\train\train.csv');
Class = categorical(labelFile.Class);
Data = [table(Class) array2table(Data_array)]; % data with labels

toc

%% Classification - 10-fold Cross Validation
% -----------------------------------------------------------
% 5 models for classification with better performances are:
%       Classification Ensmble with Subspace KNN
%       K-Nearest Neighbor
%       Cubic Kernal SVM
%       Quadratic Kernal SVM
%       Gaussian Kernal SVM
%
%   We only saved and showed the validation results for the models that
%   were relevant to our understanding and to the class (Models 1 - 5). All
%   supervised methods as our data is labeled.
%
%   Time: ~14s
% -----------------------------------------------------------

modelNames = ["Subspace KNN Ensemble","KNN","Cubic SVM","Quadratic SVM",...
    "Gaussian SVM"];

trainData = Data;       % Kfold cross validation train on entire dataset



fprintf('\n......Training classification models for 10-fold CV......\n')
tic

% Model 1: Classification Ensmble with Subspace KNN
kfoldModels{1} = fitcensemble(trainData,'Class',...
    'Method','Subspace',...     % Random subspace method
    'Learners','knn',...
    'NumLearningCycles', 20,... % Number of ensemble learning cycles
    'NPredToSample', 12);  % Number of predictors to sample for each random subspace learner

% Model 2: KNN
kfoldModels{2} = fitcknn(trainData,'Class');

% Model 3: Cubic Kernal SVM
cubicSVM = templateSVM(...
    'KernelFunction', 'polynomial', ...
    'PolynomialOrder', 3, ...   % Cubic: order = 3
    'KernelScale', 'auto', ...  % Automatically select a scale factor using a heuristic procedure
    'Standardize', true);       % Centers and scales each feature, save time for training
% Fit to multiclass SVM
kfoldModels{3} = fitcecoc(trainData,'Class','Learners',cubicSVM);   

% Model 4: Quadratic Kernal SVM
quadSVM = templateSVM(...
    'KernelFunction', 'polynomial', ...
    'PolynomialOrder', 4, ...   % Quadratic: order = 4
    'KernelScale', 'auto', ...
    'Standardize', true);
kfoldModels{4} = fitcecoc(trainData,'Class','Learners',quadSVM);

% Model 5: Gaussian Kernal SVM
guassSVM = templateSVM(...
    'KernelFunction', 'gaussian', ...
    'PolynomialOrder', [], ...
    'KernelScale', 'auto', ...
    'Standardize', true);
kfoldModels{5} = fitcecoc(trainData,'Class','Learners',guassSVM);
toc

% Predict classes for the test data set, Calculate the accuracy,
%   and plot the confusion matrices
%   Time: ~100s
fprintf('\n......Predicting and plotting confusion matrices......\n')
tic
for i=1:5
    % 10-fold cross validation
    partitionedModel = crossval(kfoldModels{i}, 'KFold', 10);
    [predictedClass, ~] = kfoldPredict(partitionedModel);
    
    % Caculate accuracy
    Accuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
    fprintf('Accuracy:  %.1f%%. CV: 10-fold. Model: %s.\n', Accuracy*100,modelNames(i))
    
    % Plot confusion matrices
    figure(i+11)
    confusionchart(Data.Class, predictedClass);
    title(strcat(modelNames(i),' Confusion Matrix 10-fold CV, Acc: ',...
        num2str(round(Accuracy,3)*100),'%'))
    saveas(figure(i+11),strcat(modelNames(i),' Confusion Matrix 10-fold CV'),'png')
end

toc
%% Classification - Train on 80% and Test on 20%
mypool = parpool();      % Set parallel computing environment

% Partition dataset into 80% trainig set and 20% test set
part = cvpartition(Data.Class, 'HoldOut',0.2);
trainData = Data(training(part),:);
testData = Data(test(part),:);


fprintf('\n......Training classification models for 80-20 Holdout CV......\n')
tic

%   Model Creation and  Holdout validation:
% ------------------------------------------------------------------------
%   Same models for classification are used as above
%   but here only trianData are used for training.

Models{1} = fitcensemble(trainData,'Class',...
    'Method','Subspace',...
    'Learners','knn',...
    'NumLearningCycles', 20,...
    'NPredToSample', 12);

Models{2} = fitcknn(trainData,'Class');

cubicSVM = templateSVM(...
    'KernelFunction', 'polynomial', ...
    'PolynomialOrder', 3, ...
    'KernelScale', 'auto', ...
    'Standardize', true);
Models{3} = fitcecoc(trainData,'Class','Learners',cubicSVM,...
    'Options',statset('UseParallel',true));

quadSVM = templateSVM(...
    'KernelFunction', 'polynomial', ...
    'PolynomialOrder', 4, ...
    'KernelScale', 'auto', ...
    'Standardize', true);
Models{4} = fitcecoc(trainData,'Class','Learners',quadSVM,...
    'Options',statset('UseParallel',true));

guassSVM = templateSVM(...
    'KernelFunction', 'gaussian', ...
    'KernelScale', 'auto', ...
    'PolynomialOrder', [], ...
    'Standardize', true);
Models{5} = fitcecoc(trainData,'Class','Learners',guassSVM,...
    'Options',statset('UseParallel',true));
toc

% Predict classes for the test data set, Calculate the accuracy,
%   and plot the confusion matrices
fprintf('\n......Predicting and plotting confusion matrices......\n')
tic
for i=1:5
    % Predict class for testData
    predictedClass = predict(Models{i},testData);
    
    % Calculate accuracy
    Accuracy = 1 - loss(Models{i},testData);
    fprintf('Accuracy:  %.1f%%. CV: 80-20 Holdout. Model: %s.\n', Accuracy*100,modelNames(i))
    
    % Plot confusion matrices and save figures
    figure(i+16)
    confusionchart(testData.Class, predictedClass);
    title(strcat(modelNames(i),' Confusion Matrix 80-20 Holdout CV, Acc: ',...
        num2str(round(Accuracy,3)*100),'%'))
    saveas(figure(i+16),strcat(modelNames(i),' Confusion Matrix 80-20 Holdout CV'),'png')
end
toc


