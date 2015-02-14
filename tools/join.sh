#!/bin/bash

cat $1 | sed -n "$(cat $2 | awk '{printf $0"p; ";}')"
