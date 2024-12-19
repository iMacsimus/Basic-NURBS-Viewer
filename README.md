# Installation guide (Linux/OSX)

## Download

### Setup repository
Clone the repository and install submodules:
```
git clone https://github.com/iMacsimus/Basic-NURBS-Viewer.git
git submodule update --init
```

### Install Intel ISPC compiler
1) Download binaries https://ispc.github.io/downloads.html  
2) Copy ispc to /usr/local/bin or any other directory from PATH  

### Install brew (OSX)

* Download brew from https://brew.sh/

### Install cmake 
#### Linux

```
sudo apt install cmake
```

#### OSX

```
brew install cmake
```

### Install dependencies
#### Linux  

```
sudo apt install libglfw3 libglfw3-dev libsdl2-dev libomp-dev
```  

#### OSX  
```
brew install libomp glfw sdl2 embree
```
To install zenity run: 

```
brew install zenity
```

## Build
Generate build files: 

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```  

Build the application binary:  

```
cmake --build build
```

## Run
Application `NURBSViewer` will appear in `build` directory  
Run `NURBSViewer` from any working directory.  
  For example, from project root directory run:  

```
./build/NURBSViewer
```

