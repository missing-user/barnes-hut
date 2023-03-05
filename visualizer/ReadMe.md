## JS Visualizer
This visualizer is cross compiled from C++ to JavaScript using Emscripten. The command to run the build is:

    sudo emmake make -j6

The error `TypeError: expected str, bytes or os.PathLike object, not list` can be safely ignored. 

The error `error: bin/data@data does not exist` on the other hand is fatal. For some reason the bin/data folder is not created. To fix this, create the folder manually and run the build again.

The bin and data folders have the linux flags: drwxr-xr-x
