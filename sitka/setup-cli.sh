#!/bin/bash

./gradlew clean
./gradlew installDist

# Fix problem arising if eclipse is used jointly
mkdir build/xtend/test
mkdir build/blang/test

