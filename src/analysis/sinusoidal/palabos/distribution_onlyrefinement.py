#! /usr/bin/python
# Automated generation of jobs

import os

configurations = [ [ 1, 0.3, 0, 0 ], [ 1, 0.3, 0, 0.1 ], [5, 0.3, 0, 0.1], [1, 0.3, 0, 0.29], [1, 0.3, 0, 0.05], [5,0.3, 0, 0.05]]

for Refinement in range(5,11):
    os.putenv("Refinement", repr(2 ** Refinement))
    for CFG_ID in range(1,6):
        os.putenv("ChannelLength", repr(configurations[CFG_ID][0]) );
        os.putenv("ChannelAmplitude", repr(configurations[CFG_ID][3]) );
        #directory = os.getenv("NATRIUM_HOME")+"/Brenner-refinement-test/ref%i-cfg%i" %(Refinement, CFG_ID)
        #os.putenv("OUTPUT_DIR", directory)
        #if not os.path.exists(directory):
        #    os.mkdir(directory)
        os.system("qsub -N 'Brenner-Ref%i-cfg%i' Brenner-refinement-test.sh" %(2**Refinement, CFG_ID))
        #print "qsub -N 'Brenner-Ref%i-cfg%i' Brenner-refinement-test.sh" %(Refinement, CFG_ID)
