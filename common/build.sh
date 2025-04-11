build="linux_build"
if [ -n "$1" ]; then
  rm -rf $build
  mkdir $build
fi
cd $build
type="Release"
if [ -n "$2" ]; then
  type="Debug"
fi
cmake .. -DCMAKE_BUILD_TYPE=$type
make -j4
