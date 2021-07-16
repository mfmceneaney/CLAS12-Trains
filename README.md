# CLAS12-Trains

Custom trains for filtering hipo files.  (Custom CLARA wagons weren't working out.)

## Setup

First build the library with:
```bash
sed -i "s;/CLAS12-Trains;/path/to/CLAS12-Trains;g" lib/src/main/java/org/jlab/trains/*.java
./gradlew build
```
Then add the following to your startup script:
```bash
pushd /path/to/CLAS12-Trains >> /dev/null
source ./env.sh
popd >> /dev/null
```
(Use `./env.csh` for c shell.)

#

Contact matthew.mceneaney@duke.edu
