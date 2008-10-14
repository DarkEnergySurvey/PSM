# Note that this file should be "sourced" rather than executed.

if [ -n "$SHELL" ]; then
   if [ `basename $SHELL` != 'bash' ]; then
      echo "Please start a bash shell before setting up the PSM"
      return 1
   fi
fi


if [ -z "$JAVA_HOME" ]; then
    echo -n "Where is the top level directory of your Java SDK distribution?  "
    read ans
    if [ -n "$ans" ]; then 
	export JAVA_HOME="$ans"
    fi
fi

if [ -z "$DESPHOTOSTDSMOD_HOME" ]; then
    echo -n "Where is the top level directory of your DESPHOTOSTDSMOD distribution?  "
    read ans
    if [ -n "$ans" ]; then
	export DESPHOTOSTDSMOD_HOME="$ans"
    fi
fi

badconfig=0

if [ ! -d "$JAVA_HOME" ]; then
    badconfig=1
    echo "Can't find JAVA SDK directory:  $JAVA_HOME"
fi

if [ ! -d "$DESPHOTOSTDSMOD_HOME" ]; then
    badconfig=1
    echo "Can't find DES GCM directory:  $DESPHOTOSTDSMOD_HOME"
fi

if [ "$badconfig" -eq 1 ]; then
   export DESPHOTOSTDSMOD_HOME=""
   export JAVA_HOME=""
   return 1
fi

export PATH=${JAVA_HOME}/bin:${DESPHOTOSTDSMOD_HOME}/bin:${PATH}

lib=$DESPHOTOSTDSMOD_HOME/lib

CLASSPATH=$lib
for jar in $lib/*.jar
do 
    CLASSPATH=${CLASSPATH}:$jar
done
CLASSPATH=${CLASSPATH}::
export CLASSPATH

echo CLASSPATH:  $CLASSPATH
