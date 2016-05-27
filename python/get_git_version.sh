#!/bin/bash

git_hash=$(git describe --long --dirty --always)
git_date=$(git show -s --format=%ci)

cat <<EOF > version_tag.py
COMMIT = "$git_hash"
DATE = "$git_date"
EOF

