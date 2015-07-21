#!/usr/bin/python
#BSUB -W 24000
#BSUB -n 1
#BSUB -q he2
#BSUB -o log1.out
#BSUB -e log1.err

import os,subprocess

os.chdir("/home/actodd/ROMS-utils/USeast-age/age/averaging/age_above_pcline.py")
proc=subprocess.Popen(['/home/actodd/ROMS-utils/USeast-age/age/averaging/age_above_pcline.py'],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
(out,err)=proc.communicate()
