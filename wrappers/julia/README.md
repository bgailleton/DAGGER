# Julia Bindings

This folder contains an example of `DAGGER` bindings for `Julia`. This has only been tested on `linux` so far and remain quite experimental.

## installation

You will need the following dependencies:

- `cmake` (2.8 or later)
- `julia` (tested on `1.8` but should work for any version `>1`)
- `CxxWrap` package (installable within `julia`)
- a `c++` compiler compatible with c++17 (tested on `g++11`)

As a preliminary step, you need to detect a CxxWrap and Julia's path (I won't pretend I understand the details, [see here](https://github.com/JuliaInterop/CxxWrap.jl) if interested). To find that path: run the following commands:

```julia
julia
julia> using CxxWrap
julia> CxxWrap.prefix_path()
```

And save the result (let's call it PathCxxJulia for the rest of this README).

Now you are ready to compile the code (Warning: you'll have to adapt `CMakeLists.txt` if you change any file name or have a different folder structure than here):

```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=PathCxxJulia ../
cmake --build . --config Release
```

Where, as a reminder, `PathCxxJulia` is the path obtained above (I insist, because each and every time I do something like this, I get emails complaining that "PathCxxJulia is not found" by the compiler.
