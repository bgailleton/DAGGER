rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/home/bgailleton/mambaforge/envs/dagger/share/julia/artifacts/e74d2fd61bdd4d0ff365351ee632d660b2d43ec5 ../
cmake --build . --config Release
