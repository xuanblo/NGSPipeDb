import os
import sys
import datetime
import pypiper
from optparse import OptionParser

usage = '''

    Example:

        gff_stat_vistualizition.py -i test.gtf -o output_folder

'''

parser = OptionParser(usage)
parser.add_option('-n', '--pipename', dest='pipename', help='name for piper', action = "store", type="string")
parser.add_option("-o","--output-dir", dest="outfolder", help="save file to this directory", action = "store", type="string")
parser.add_option("-c","--cpu", dest="cpu", help="save file to this directory", action = "store", type="string")
parser.add_option("-m","--memmory", dest="mem", help="how many memmory to use", action = "store", type="string")
parser.add_option("-t","--log", dest="trinitylog", help="log file", action = "store", type="string")
parser.add_option("-s","--sample", dest="sam", help="trinity sample file", action = "store", type="string")


(options, args)=parser.parse_args()

if os.path.exists(options.outfolder):
    os.removedirs(options.outfolder)

pm = pypiper.PipelineManager(name=options.pipename, outfolder=options.outfolder)

pm.timestamp("start assembly using trinity!")

command = "Trinity --seqType fq --max_memory {} --output {} --samples_file {} --CPU {}".format(options.mem, options.outfolder, options.sam, options.cpu)

target_file = options.trinitylog

pm.run(command, target_file)

pm.stop_pipeline()
