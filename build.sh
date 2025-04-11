build="linux_build"
comm="common"
if [ -n "$1" ]; then
  cd $comm
  sh build.sh 1 $2
  cd ../
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
