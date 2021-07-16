#!/bin/bash
export JM_OPTS="-Xms32m -Xmx1500m -XX:+UseG1GC"
java -cp $CLASSPATH\:lib/build/libs/lib.jar org.jlab.trains.Train $@
