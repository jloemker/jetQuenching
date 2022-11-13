#!/bin/sh

#copy file in tutorial section
cp correlationV0jet.cxx /home/johannalomker/alice/O2Physics/Tutorials/src/correlationV0jet.cxx 
#enter installation directory 
cd /home/johannalomker/alice/sw/BUILD/O2Physics-latest/O2Physics/
pwd
#rebuild whole tutorial section via ninja install
eval $(alienv printenv O2Physics/latest ninja/latest)
ninja install Tutorials/all 

