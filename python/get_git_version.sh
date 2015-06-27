#!/bin/bash

git_hash=$(git rev-parse HEAD)

cat <<EOF > version_tag.py
def print_version():
    print "\n VERSION INFO"
    print " Git commit: $git_hash"
EOF

