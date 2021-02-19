#!/bin/bash

echo "Home dir: $HOME"
echo "Creating temp dir"
tmp_dir=$(mktemp -d -t test-XXXXXXX)
echo "Created temp dir $tmp_dir"
echo "Removing $tmp_dir"
rm -rf "$tmp_dir"

exit 0
