set CORRLIBPATH=`correction config --libdir`
echo $CORRLIBPATH
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$CORRLIBPATH
