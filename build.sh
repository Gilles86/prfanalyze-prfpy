#! /bin/bash

[ -z "$USERTAG" ] && USERTAG=niklasmueller
[ -z "$SOLVER"  ] && SOLVER=prfpy
[ -z "$VERSION" ] && VERSION=latest

SCRIPTPATH="$( cd "$(dirname "$0")" && pwd -P )"

echo "Building $USERTAG/analyze-$SOLVER:$VERSION"
echo "  Path: $SCRIPTPATH"
echo "------------------------------------------------------------"
echo ""

docker build "$@" --tag "$USERTAG/analyze-$SOLVER:$VERSION" "$SCRIPTPATH"