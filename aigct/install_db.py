import argparse
from aigct.install_util import init_db


def main():
    parser = argparse.ArgumentParser(description="Install Database")
    parser.add_argument('--confdir', type=str,
                        help='Config file directory')
    args = parser.parse_args()

    if not args.confdir:
        raise Exception(
            "Must specify --confdir argument")

    init_db(args.confdir)


if __name__ == "__main__":
    main()
