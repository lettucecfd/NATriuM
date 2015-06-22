#! /usr/bin/python
# Automated generation of jobs

import os

os.putenv("Ma", repr(0.01))
os.putenv("Gamma", repr(1))
os.putenv("CURRENT_DIR", os.path.dirname(os.path.realpath(__file__)))

for orderFE in range(2,10,2):
    os.putenv("orderFE", repr(orderFE))
    for Refinement in range(2,5):
        os.putenv("Refinement", repr(Refinement))
        for CFG_ID in range(1,6):
            os.putenv("CFG_ID", repr(CFG_ID))
	    directory = os.getenv("NATRIUM_HOME")+"/Brenner-refinement-test/ref%i-p%i-cfg%i" %(Refinement, orderFE, CFG_ID)
            os.putenv("OUTPUT_DIR", directory)
            if not os.path.exists(directory):
                os.mkdir(directory)
            os.system("qsub -N 'Brenner-Ref%i-p%i-cfg%i' Brenner-refinement-test.sh" %(Refinement, orderFE, CFG_ID))
