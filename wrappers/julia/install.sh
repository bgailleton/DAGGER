rm -rf build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/home/bgailleton/mambaforge/envs/dagger/share/julia/artifacts/dadee2d976c39012e3417bdcb3a73e45adb7ae7c ../
cmake --build . --config Release
cd ../
cp build/lib/libjudagger* demo/
