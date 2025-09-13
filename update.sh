SCRIPT_PATH="$(readlink -f "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
cd $SCRIPT_DIR/MedialAxisCGAL/build
cmake ../src/ -DIPE_PREFIX=/usr/local/ipe
cmake --build .
sudo mv lib/medialaxis_ipelet.so /usr/local/lib/ipe/7.3.1/ipelets/medialaxis.so
cd -
