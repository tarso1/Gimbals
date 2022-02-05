# Gimbals

## Table of contents
* [General info](#general-info)
* [Technologies](#technologies)
* [Setup](#setup)

## General Info
 This program simulates rotation movements through a set of three gimbals connected, the inner, the middle, and the outer gimbal. Gimbals are like frames that can rotate independently. The inner gimbal may represent any object rotating in three dimensions, such as a car, a rocket, an airplane, a satellite, or even a hand. The composition of three independent rotations can reach any orientation in space. To achieve this, the program has been made using quaternions to represent the orientation and the quaternions kinematics to produce the movement. At the end of each simulation, the program shows the final quaternion, i.e., the orientation of the inner gimbal with respect to an inertial reference frame, as well the Euler angles associated using the convention Yaw, Pitch, and Roll. 
 ## Technologies
 The technologies used to build this program were:
 * Python 3.9
 * Numpy
 * Scipy
 * Matplotlib
 * Tkinter
 
 ## Setup
 To run this project you need to execute the python script named "program_gimbals.py". For this you need to have installed Python 3.9 or more recent versions. To check your python version open a linux terminal or the cmd for Windows and type:
 ```
 python3 -V
 ```
 To install pip:
 ```
 sudo apt install pip3
```
 And the libraries numpy, scipy, matplotlib:
 ```
 pip3 install matplotlib
 pip3 install numpy
 pip3 install scipy
 ```
 To run the program, go to the folder where is the script "program_gimbals" and execute:
 ```
 python3 program_gimbals.py
 ```
![Algorithm schema](./images/window.png)
