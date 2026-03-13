import grunner as grun
import pprint

pp = pprint.PrettyPrinter()


# Modify the defaults
grun.data["finaltime"] = 10.

# Print out the modifications 
pp.pprint(grun.data)

grun.data["NX"] =32
grun.data["outputfiletag"]=grun.getdefault_filename("grunexample")
grun.run(dry_run=False, ncpus="1")
