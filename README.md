
# Aparapi accelerated Java Barnes Hut code

A simple Java code illustrating use of Barnes Hut algorithm for N-body simulation of a group of unit mass bodies - "stars" - interacting through the force of gravity.

This is a parallel version designed to run on systems with OpenCL-compatible graphics acceleration.  It does this through the Java API [Aparapi](www.aparapi.com).

## Getting Started

### Prequisites

Oracle Java 8 JDK and Maven 3 for the build (note Oracle JDK is specified for Aparapi).  A computer with OpenCL available (often available by default with the drivers of the graphics card).

### Installing and running

Download distribution from GitHub.

In the project folder, build by:
```
  $ mvn package
```
Run by:
```
  $ java -jar target/aparapi-test-1.0-SNAPSHOT-jar-with-dependencies.jar
```
A Java graphics window should appear to display current state of
simulation.  Monitoring output including profiling information will be printed at the terminal.

Early in the monitoring output of the program, you should see a message about "Device Usage by Kernel".  If you are successfully running on the graphics card this may be followed by, e.g.:
```
org.hpjava.KernelTree:
        using Intel<GPU>
```
Any more complicated message here suggests the GPU may not have been successfully invoked - Aparapi typically then falls back to a Java Thread Pool.

The simulation will continue running until the graphics window is closed or the program is killed at the terminal.

## Disclaimer

Although the logic in this code is believed to be a correct implementation of Barnes-Hut, parameters including the time step and opening angle have not been tuned to guarantee accuracy of the simulation.

