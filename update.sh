SCRIPT_PATH="$(readlink -f "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
cd $SCRIPT_DIR/MedialAxisCGAL/build
cmake ../src/ -DIPE_PREFIX=/usr/local/ipe -DBUILD_CLI=false
cmake --build .
sudo cp $SCRIPT_DIR/MedialAxisCGAL/build/medialaxis.so /usr/local/lib/ipe/7.3.1/ipelets/medialaxis.so
sudo cp $SCRIPT_DIR/MedialAxisLua/medialaxis.lua /usr/local/lib/ipe/7.3.1/ipelets/.
cd -
