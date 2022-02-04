#! /bin/bash

[ -z "$USERTAG" ] && USERTAG=niklasmueller01
[ -z "$SOLVER"  ] && SOLVER=prfpy
[ -z "$VERSION" ] && VERSION=latest

SCRIPTPATH="$( cd "$(dirname "$0")" && pwd -P )"

echo "Building $USERTAG/prfanalyze-$SOLVER:$VERSION"
echo "  Path: $SCRIPTPATH"
echo "------------------------------------------------------------"
echo ""

docker build "$@" --tag "$USERTAG/prfanalyze-$SOLVER:$VERSION" "$SCRIPTPATH" --network=host