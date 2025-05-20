import argparse
from aigct.install_util import check_install


def main():
    parser = argparse.ArgumentParser(description="Check Install")
    parser.add_argument('--confdir', type=str,
                        help='Config file directory')
    args = parser.parse_args()

    if not args.confdir:
        raise Exception(
            "Must specify --confdir argument")

    check_install(args.confdir)


if __name__ == "__main__":
    main()
