#! /bin/csh -f
# Note that this file should be "sourced" rather than executed.

if (! $?JAVA_HOME) then
    echo -n "Where is the top level directory of your Java SDK distribution?"
    set ans = "$<"
    if ($ans != "") setenv JAVA_HOME $ans
endif

if (! $?DESPHOTOSTDSMOD_HOME) then
    echo -n "Where is the top level directory of your DESPHOTOSTDSMOD distribution?"
    set ans = "$<"
    if ($ans != "") setenv DESPHOTOSTDSMOD_HOME $ans
endif

set installation = ("DESPhotoStdsMod" "Java SDK")
set installdir = ($DESPHOTOSTDSMOD_HOME $JAVA_HOME)
unset badconfig
foreach i (1 2)
   if (! -d $installdir[$i]) then
       set badconfig = 1
       echo "Can't find $installation[$i] directory: $installdir[$i]"
   endif
end
if ($?badconfig) then
   unsetenv DESPHOTOSTDSMOD_HOME
   unsetenv JAVA_HOME
   exit(1)
endif
echo DESPHOTOSTDSMOD_HOME:  $DESPHOTOSTDSMOD_HOME
echo JAVA_HOME:  $JAVA_HOME

set path = ($JAVA_HOME/bin $DESPHOTOSTDSMOD_HOME/bin $path)

set lib = $DESPHOTOSTDSMOD_HOME/lib

setenv CLASSPATH .
setenv CLASSPATH ${CLASSPATH}:$lib
foreach jar ($lib/*.jar)
    setenv CLASSPATH ${CLASSPATH}:$jar
end

echo CLASSPATH:  $CLASSPATH
