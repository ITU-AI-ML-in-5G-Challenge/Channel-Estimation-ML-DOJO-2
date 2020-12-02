# Channel-Estimation-ML-DOJO-2


# The problem statement can be found here: 
[Challenge: Site-specific channel estimation with hybrid MIMO architectures](https://research.ece.ncsu.edu/ai5gchallenge/) 
The train and evaluation data can also be downloaded from that link.
# Requirements:
  The code has been developed using Matlab 2019b.
# Code:
* TD1SNR1_example_for_test_data.m is an example of how to obtain the predictions for test data for the test dataset 1- SNR1.
* ReconstructChannel.m is the function with our algorithm to obtain the channels.
   Whitening.m is a subfunction used to whitten the noise on ReconstructChannel.m.
   
* Train.m is the script used for training to obtain the detection threshold.
   ReconstructChannel_train.m is the function similar to ReconstructChannel.m used during the training script.

# Results:
The test channels can be downloaded from [this drive folder] (https://1drv.ms/u/s!AouYq8xpaglhhZgVADXpc5YFDm2G3w?e=cVmsGW)
