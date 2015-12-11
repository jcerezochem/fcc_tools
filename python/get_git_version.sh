#!/bin/bash

git_hash=$(git describe --long --dirty --always)
git_date=$(git show -s --format=%ci)

cat <<EOF > version_tag.py
def print_version():
    print "\n VERSION INFO"
    print " Git commit : $git_hash"
    print " Commit date: $git_date"
EOF

