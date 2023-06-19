# Hybrid NEMO proofs of concept - Oriol Duran

This repository contains the code generated during the development of my TFG research project, which consisted of the implementation of three proofs of concept of the NEMO model with hybrid parallelism at the bottleneck.

The different parallel strategies correspond to the branches as follows:

- loops strategy: loops
- paralllel tracer loop strategy: Parallel_tradvfct_comms
- nested stratedy: Reverse_No_Collapse

## Compiling the proofs of concept
This repository only contains the code that needed to be modified in order to implement said parallelism.

To compile and run the model for testing you can follow the guide at https://sites.nemo-ocean.io/user-guide/install.html#download-and-install-the-nemo-code to download the base NEMO model (note that this project was developed on the NEMO 4.2.0 version).

After downloading the base NEMO model code:

1. copy this repository in the same directory where you downloaded the base model and let the new files overwrite the ones from the base model.

2. proceed with the instructions at https://sites.nemo-ocean.io/user-guide/install.html#download-and-install-the-nemo-code to compile the code, be sure to add the OpenMP flag to your arch file before executing the compilation command (-qopenmp in our archfiles, but may vary depending on which OpenMP implementation you are using)

## Running tests with the proofs of concept
Once the compilation is completed, the steps to prepare and execute any test are the same as they would be for the base NEMO (https://sites.nemo-ocean.io/user-guide/install.html#preparing-an-experiment)
