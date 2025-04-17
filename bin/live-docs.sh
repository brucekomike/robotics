#!/usr/bin/env bash

if (! command -v sphinx-autobuild ) ; then
  pip install -r ./docs/requirements.txt
fi
sphinx-autobuild ./docs/source ./docs/build