#!/bin/bash
export JM_OPTS="-Xms32m -Xmx1500m -XX:+UseG1GC"
java -cp $CLASSPATH\:$C12TRAINS/lib/build/libs/lib.jar $JAVA_OPTS $JM_OPTS org.jlab.trains.LambdaTrain $@
