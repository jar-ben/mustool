#!/bin/bash
#
# File:  update_defs.sh
# Author:  mikolas
# Created on:  07 Apr 2017 18:20:12
# Copyright (C) 2017, Mikolas Janota
#
TF=/tmp/defs_${RANDOM}
GITHEAD=`git rev-parse HEAD`
echo '#define GITHEAD "'${GITHEAD}'"' >${TF}
if [ -e defs.h ]; then
    if diff -q defs.h ${TF} >/dev/null; then
        echo defs remain;
    else
        /bin/cp ${TF} defs.h
    fi
else
    cp ${TF} defs.h
fi
rm ${TF}
