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
