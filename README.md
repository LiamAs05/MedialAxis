# Medial Axis Ipelet
An [ipelet](https://ipe.otfried.org/) that computes the medial axis of a selected polygon.

![Example](.rsrc/example.png?raw=true "Medial Axis of a Simple Polygon")

## Compiling (Linux Instructions)

Make sure you have `ipe` installed with its [dependencies](https://github.com/otfried/ipe/blob/master/doc/install.txt).

Make sure you have `cmake` (`sudo apt install cmake`)

After cloning the project, just run `./update.sh` from the projects root directory, it should automatically place the ipelet at the first `ipelets` directory it manages to find.

NOTE: if you want to specify the ipelets dir manually, just run `./update.sh <directory>`.
