from aigct.init_app import main
import sys

sys.argv = ["junk1.py", "--confdir", "/tmp/conf", "--dbdir", "/tmp/db", "--logdir", "/tmp/log", "--outdir", "/tmp/out"]
main()