import urllib2
import json
import sys
import os
import gzip 
from . import __file__

def download(args):
    """download Pfam HMMs and write download destination into config file"""
    attempts = 0
    if not args.local: 
        while attempts < 3:
            try:
                if not os.path.exists(args.download):
                    sys.stderr.write("directory %s doesn't exists; please create first\n" % args.download)
                    sys.exit(0)
                response = urllib2.urlopen("ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz", timeout = 5)
                with open(os.path.join(args.download, "Pfam-A.hmm.gz"), 'w' ) as f:
                    CHUNK = 1000000
                    while True:
                        chunk = response.read(CHUNK)
                        if not chunk:
                            break
                        else:
                            f.write(chunk)
                with gzip.open(os.path.join(args.download, "Pfam-A.hmm.gz"), 'rb') as zf:
                    with open(os.path.join(args.download, "Pfam-A.hmm"), 'wb') as out_f:
                        for l in zf:
                            out_f.write(l)
                break
            except urllib2.URLError as e:
                attempts += 1
                print e
    with open(os.path.abspath(os.path.dirname(__file__)) + "/" + "config.json", 'w') as config:
        if not os.path.exists(os.path.join(args.download, "Pfam-A.hmm")):
            sys.exit("something went wrong; make sure %s contains Pfam-A.hmm" % args.download)
        config.write(json.dumps({"hmms": os.path.abspath(os.path.join(args.download))}))
