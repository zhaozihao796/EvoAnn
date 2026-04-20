import logging

from .cli import create_parser
import sys


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.StreamHandler()]
    )
    parser = create_parser()
    # if len(sys.argv) == 1:
    #     parser.print_help()
    #     sys.exit(0)
    if len(sys.argv) == 1 or '-h' in sys.argv or '--help' in sys.argv:
        if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1] in ['-h', '--help']):
            parser.print_help()
            sys.exit(0)
    try:
        args = parser.parse_args()
        if hasattr(args, 'func'):
            args.func(args)
        else:
            parser.print_help()
    except KeyboardInterrupt:
        sys.exit(130)
    except Exception as e:
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
