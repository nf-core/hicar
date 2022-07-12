#!/usr/bin/env python

#########################################
# Author: Jianhong Ou
#########################################

import Bio.Restriction as biorst
import argparse
import sys

def parse_args(args=None):
    Description = 'Get restriction enzyme cutting sites'
    Epilog = """Example usage: python restriction_enzyme_cutsite.py <enzyme>"""
    argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    ## REQUIRED PARAMETERS
    argParser.add_argument('ENZYME', help="enzyme name")
    return argParser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    if args.ENZYME.lower()=='mnase':
        sys.stdout.write('mnase 0')
    else:
        try:
            en = getattr(biorst, args.ENZYME)
        except AttributeError:
            raise ValueError("Unknown enzyme name: {}".format(args.ENZYME))

        sys.stdout.write(en.site+' '+str(en.fst5))

    sys.stdout.flush()
    return(0)


if __name__ == "__main__":
    sys.exit(main())
