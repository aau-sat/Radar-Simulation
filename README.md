# Radar-Simulation

Code used for the 3D driving scenario in the parking garage simulation, as described in the paper:

A holistic, high-fidelity and noise aware FMCW radar simulation framework: from waveform to point cloud

submitted to the IEEE Sensor Journal.

The main code is "RadarSimulation_3D_driving". The code connects to Unity via TCP/IP and provides:
% - time domain data with noise
% - range FFT
% - doppler FFT
% - Range Doppler Map
% - Range Angle Map
% - 3D radar point cloud

Phase Noise data are not given, by request of the provider. As such, the function that generates DPN samples is also not provided in this repository.

