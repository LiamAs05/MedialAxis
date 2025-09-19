#!/usr/bin/env bash
set -e

SCRIPT_PATH="$(readlink -f "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"

# If a path is given, use it; otherwise search
if [ -n "${1:-}" ]; then
    IPELETS_DIR="$1"
else
    IPELETS_DIR="$(find / -name ipelets 2>/dev/null | head -n 1)"
fi

echo "Using ipelets directory: $IPELETS_DIR"
IPELET_FILES="$SCRIPT_DIR/MedialAxisCGAL/build/medialaxis.so $SCRIPT_DIR/MedialAxisLua/medialaxis.lua"

mkdir -p "$SCRIPT_DIR/MedialAxisCGAL/build"
cd "$SCRIPT_DIR/MedialAxisCGAL/build"
cmake ../src/ -DIPE_PREFIX=/usr/local/ipe -DBUILD_CLI=false
cmake --build .
sudo cp $IPELET_FILES "$IPELETS_DIR"
echo -n "Changing directory back to: "
cd -
