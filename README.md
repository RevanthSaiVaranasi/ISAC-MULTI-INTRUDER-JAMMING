# ISAC Multi-Intruder Detection and Jamming Framework

This repository contains a MATLAB implementation of a sensing-aware jamming framework for multi-intruder detection and interception in Integrated Sensing and Communication (ISAC) systems.

The project is divided into three major algorithms:

 Stage-1: Multi-intruder detection using Capon beamforming and CRB-based uncertainty interval generation.
 Stage-2: Robust beamforming optimization to satisfy worst-case Jamming-to-Noise Ratio (JNR) constraints under transmit power limitations.
 Stage-3: Intruder scheduling scheme that identifies feasible jamming groups and sequentially intercepts all detected intruders when simultaneous jamming is not possible.
 
Features

- Multi-target direction-of-arrival estimation
- CRB-based uncertainty modeling
- Robust beamforming design
- JNR-constrained jamming optimization
- Multi-round intruder scheduling
- MATLAB-based simulation framework

Repository Structure

- "Stage_1.m" – Detection and angle estimation
- "Stage_2.m" – Robust beamforming optimization
- "Stage_3.m" – Intruder scheduling
- "steering_vector.m" – Utility function for antenna array steering vectors

Author

Revanth Sai Varanasi
B.Tech, Electronics and Communication Engineering
National Institute of Technology Warangal
