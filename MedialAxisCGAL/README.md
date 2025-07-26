# Compiling The Project

You should have a `src` directory containing the `MedialAxis.cpp` and `CMakeLists.txt` files only.

To configure a VS Project, you can follow the steps below, make sure `cmake` is installed.

## Installing CGAL On Windows

### Install vcpkg

```cmd
cd <CHOOSE_VCPKG_INSTALL_PATH>
git clone https://github.com/microsoft/vcpkg
cd vcpkg
.\bootstrap-vcpkg.bat
```

### Set VCPKG_ROOT

```cmd
setx VCPKG_ROOT "<VCPKG_INSTALL_PATH>\vcpkg"
```

### Install CGAL

This takes way too long...

```cmd
vcpkg.exe install yasm-tool:x86-windows
vcpkg.exe install cgal
```

## Configuring The Project

```cmd
cd <DIR_WITH_MEDIAL_AXIS_CPP_FILE>
mkdir build
cd build
cmake-gui ..
```

- Specify the source code path (`<DIR_WITH_MEDIAL_AXIS_CPP_FILE>`).
- Specify `Visual Studio 17 2022` as the `Generator`.
- Specify `x64` as the `Optional Platform`.
- Select `Specify toolchain file for cross compilation` and choose the `vcpkg.cmake` file, located in `<VCPKG_INSTALL_PATH>/scripts/buildsystems/vcpkg.cmake`

