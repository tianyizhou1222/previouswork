%% Exported from Jupyter Notebook
% Run each section by placing your cursor in it and pressing Ctrl+Enter
disp('Launching...')

%% Code Cell[1]:

addpath(".")

%% Code Cell[2]:

% Set this to the directory where the data is stored
dataFolder = "~/Music/urbansound8k";

%% Code Cell[3]: Get labels for data set

ads = audioDatastore(strcat(dataFolder,"/Train"));
[~, names, ~] = cellfun(@fileparts, ads.Files, 'UniformOutput', false);
names = cellfun(@str2num, names);
[~, idxs] = sort(names);
ads.Files = ads.Files(idxs);
labels = readtable(strcat(dataFolder,"/train.csv"));
ads.Labels = categorical(labels.Class);

%% Code Cell[4]:

[sTrain, sVal] = splitEachLabel(ads, 0.8);

% # Neural Networks
% Create spectrograms from the data. The dimension of the output data will then be 40x396. Shorter audio clips are padded equally on both sides, see `spectrograms.m` for details.

%% Code Cell[20]:

segmentDuration = 4;
frameDuration = 0.05;
hopDuration = 0.010;
numBands = 40;
reset(sTrain);
reset(sVal);
disp('Generating Spectrograms...')
epsilon = 1e-6; % Added so that log doesn't encounter 0
XTrain = log10(spectrograms(sTrain, segmentDuration, frameDuration, hopDuration, numBands)+epsilon);
XVal = log10(spectrograms(sVal, segmentDuration, frameDuration, hopDuration, numBands)+epsilon);


% Create categorical vectors for labels
YTrain = categorical(sTrain.Labels);
YVal = categorical(sVal.Labels);


% Visualize some training data
n = randi([0,500], 1, 1);
figure
pcolor(XTrain(:,:,:,n))
title(strrep(char(YTrain(n,:)), '_', ' '))
shading flat
[samps, sampfreq] = audioread(sTrain.Files{n});
sound(samps, sampfreq)

% ## Conv Net from Scratch
% This section uses a simple convolutional neural network to classify spectrograms and achieves a validation set accuracy of 95.86% using holdout validation set of 20% of the data. The network for this section is largely the same as one from a speech recognition tutorial by Mathworks [here](https://www.mathworks.com/help/deeplearning/examples/deep-learning-speech-recognition.html?s_tid=mwa_osa_a). This CNN is trained using the [Adam](https://arxiv.org/abs/1412.6980) optimizer, and begins to overfit starting at about 7 epochs using the learning rate $3*10^{-3}$.

sz = size(XTrain);
specSize = sz(1:2);
imageSize = [specSize 1];
numClasses = numel(categories(YTrain));

timePoolSize = ceil(imageSize(2)/8);
dropoutProb = 0.2;
numF = 12;
layers = [
    imageInputLayer(imageSize, 'Name', 'Input_Layer')

    convolution2dLayer(3,numF,'Padding','same', 'Name', 'Conv_1')
    batchNormalizationLayer('Name', 'BN_1')
    reluLayer('Name', 'Relu_1')

    maxPooling2dLayer(3,'Stride',2,'Padding','same', 'Name', 'MaxPool_1')

    convolution2dLayer(3,2*numF,'Padding','same', 'Name', 'Conv_2')
    batchNormalizationLayer('Name', 'BN_2')
    reluLayer('Name', 'Relu_2')

    maxPooling2dLayer(3,'Stride',2,'Padding','same', 'Name', 'MaxPool_2')

    convolution2dLayer(3,4*numF,'Padding','same', 'Name', 'Conv_3')
    batchNormalizationLayer('Name', 'BN_3')
    reluLayer('Name', 'Relu_3')

    maxPooling2dLayer(3,'Stride',2,'Padding','same', 'Name', 'MaxPool_3')

    convolution2dLayer(3,4*numF,'Padding','same', 'Name', 'Conv_4')
    batchNormalizationLayer('Name', 'BN_4')
    reluLayer('Name', 'Relu_4')
    convolution2dLayer(3,4*numF,'Padding','same', 'Name', 'Conv_5')
    batchNormalizationLayer('Name', 'BN_5')
    reluLayer('Name', 'Relu_5')

    maxPooling2dLayer([1 timePoolSize], 'Name', 'MaxPool_4')

    dropoutLayer(dropoutProb, 'Name', 'Dropout_1')
    fullyConnectedLayer(numClasses, 'Name', 'FC_1')
    softmaxLayer('Name', 'Softmax_1')
    classificationLayer('Name', 'Classification')];


% Plot the layer graph of the network
plot(layerGraph(layers))

% Training options
miniBatchSize = 128;
valFrequency = floor(numel(YTrain)/miniBatchSize);
options = trainingOptions('adam',...
    'InitialLearnRate',3e-3, ...
    'MaxEpochs',25, ...
    'MiniBatchSize',miniBatchSize, ...
    'Shuffle','every-epoch', ...
    'Plots','training-progress', ...
    'Verbose',false, ...
    'ValidationData',{XVal,YVal}, ...
    'ValidationFrequency',valFrequency, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.1, ...
    'LearnRateDropPeriod',20, ...
    'ExecutionEnvironment', 'gpu');

% Train network on the GPU!
disp('Training Feedforward Network...')
trainedNet = trainNetwork(XTrain, YTrain, layers, options);

% Get training and validation loss from trained model
YValPred = classify(trainedNet, XVal);
valError = mean(YValPred ~= YVal);
YTrainPred = classify(trainedNet, XTrain);
trainError = mean(YTrainPred ~= YTrain);
disp("Training Error: "+ trainError*100+"%")
disp("Validation Error: "+ valError*100+"%")

% # confusion matrix
figure('Units','normalized','Position',[0.2 0.2 0.5 0.5]);
cm = confusionchart(YVal,YValPred);
cm.Title = 'Confusion Matrix for Validation Data';


% ## Visualize Early Conv Layer features - nice and simple shapes
chans = 1:12
I = deepDreamImage(trainedNet, 2, chans, 'PyramidLevel', 1);
figure
I = imtile(I,'ThumbnailSize',[64 64]);
imshow(I)
title(['Visualization of Layer conv_1 Features'])

% ## Visualize Late Conv Layer - very complicated shapes! textures really.
chans = 1:48
I = deepDreamImage(trainedNet, 17, chans, 'PyramidLevel', 1);
figure
I = imtile(I,'ThumbnailSize',[64 64]);
imshow(I)
title(['Visualization of Layer conv_5 Features'])

