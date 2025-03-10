import argparse
from aigct.install_util import init_config_file, init_db


def main():
    parser = argparse.ArgumentParser(description="Init App")
    parser.add_argument('--confdir', type=str,
                        help='Config file directory')
    parser.add_argument('--dbdir', type=str,
                        help='Directory where db files are to be installed')
    parser.add_argument('--logdir', type=str,
                        help='Directory where log files are to be written')
    parser.add_argument(
        '--outdir', type=str,
        help='Directory where analysis output files are to be written')
    args = parser.parse_args()

    if (not args.confdir or not args.dbdir or not args.logdir or
            not args.outdir):
        raise Exception(
            "Must specify --confdir, --dbdir, --logdir, and --outdir arguments")

    init_config_file(args.confdir, args.dbdir, args.outdir,
                     args.logdir)
    # init_db(args.confdir)


if __name__ == "__main__":
    main()
