%% Exported from Jupyter Notebook
%
% Authors: Zachary Neveu | Tanyi Zhou | Christian Grenier
% NOTE: Run feedforward_nn.m first to generate spectrograms
%
%% Markdown Cell:
% # Pre-Trained Resnet
% Residual networks (Resnets) consist of a series of "residual blocks" of layers, and bypass connections to skip these layers. The idea of these types of networks is that each residual block can correct for the error in the calculation of the previous res block, allowing the network to operate on errors rather than full results. Resnets generally show improved performance over traditional convolutional networks.
% 
% This section also utilizes transfer learning, making use of a model that has already been trained to recognize images from the Image Net dataset. While spectrograms are very different from real-world images, the early layers of the network should be similar, capturing low-level features such as edges. In order to re-train part of this network, the final layer is replaced to get a 10 class classifier for our dataset. The network is then trained as if the weights had been randomly initialized. One way this process could be improved is to set a variable learning rate so that early layers in the network adjust slower than the final layers of the network whose weights are largely useless for our purposes.
% 
% In order to use the pre-trained resnets that Mathworks provides, input images must be 224x224 pixels. We must then resize our spectrograms to be fed into this network. The network also expects color images. We add color to our spectrogram representations using the `pcolor()` function. Research has shown (CITE THIS) that adding color can actually improve network performance. Reprocess spectrograms as images with same size as pretrained resnet (WHYYYYY doesn't matlab handle varying input sizes for you???)
% 
% This method achieves 95.86% accuracy on the same validation set as the previous conv net. This is very similar in terms of performance, however the pre-trained resnet is able to reach this accuracy in only 3 epochs.

%% Markdown Cell:
% ![performance](pres/resnet_3_epochs.png)

%% Code Cell[22]:

% (Only if data pre-generated)
%load("resnet") 

%% Code Cell[23]:

disp('Resizing Validation spectrograms, be patient...')
XValIms = specs2Ims(XVal, [224 224]);
disp('Resizing Training spectrograms, be extra patient...')
XTrainIms = specs2Ims(XTrain, [224 224]);

%% Code Cell[24]:

%save("resnet", '-v7.3')

%% Code Cell[55]:

net = resnet18;
inputSize = net.Layers(1).InputSize

%% Markdown Cell:
% Use Matlab [example](https://www.mathworks.com/help/deeplearning/examples/train-deep-learning-network-to-classify-new-images.html) supporting function to get last two layers of network

%% Code Cell[56]:

lgraph = layerGraph(net);
[learnableLayer,classLayer] = findLayersToReplace(lgraph);

%% Code Cell[57]:

numClasses = 10;
newLearnableLayer = fullyConnectedLayer(numClasses, ...
    'Name','new_fc', ...
    'WeightLearnRateFactor',10, ...
    'BiasLearnRateFactor',10);

%% Markdown Cell:
% Replace the last pre-trained layer and classification with untrained layers based on the number of classes we require

%% Code Cell[58]:

lgraph = replaceLayer(lgraph,learnableLayer.Name,newLearnableLayer);
newClassLayer = classificationLayer('Name','new_classoutput');
lgraph = replaceLayer(lgraph,classLayer.Name,newClassLayer);

%% Markdown Cell:
% Optionally, freeze initial layers of resnet here (kind of complicated compared to fastai)

%% Markdown Cell:
% Use same optimizer as before to train

%% Code Cell[59]:

miniBatchSize = 128;
valFrequency = floor(numel(YTrain)/miniBatchSize);
options = trainingOptions('adam',...
    'InitialLearnRate',3e-4, ...
    'MaxEpochs',3, ...
    'MiniBatchSize',miniBatchSize, ...
    'Shuffle','every-epoch', ...
    'Plots','training-progress', ...
    'Verbose',false, ...
    'ValidationData',{XValIms,YVal}, ...
    'ValidationFrequency',valFrequency, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.1, ...
    'LearnRateDropPeriod',2, ...
    'ExecutionEnvironment', 'gpu');

%% Code Cell[60]:

disp('Training Resnet...')
trainNetwork(XTrainIms, YTrain, lgraph, options)

%% Markdown Cell:
% Fin.
