#!/usr/bin/env python3
"""Compile the manuscript with Tectonic, keeping temporary files out of the repository."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
PDF_NAME = 'generalized-multiplication-tables.pdf'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='tectonic',
                        help='Tectonic executable name or path (default: tectonic)')
    parser.add_argument('--offline', action='store_true',
                        help='use only already cached standard TeX resources')
    args = parser.parse_args()
    compiler = shutil.which(args.compiler)
    if compiler is None:
        parser.error('Tectonic was not found. Install it or pass --compiler /path/to/tectonic.')
    compiler = str(Path(compiler).resolve())

    with tempfile.TemporaryDirectory(prefix='multiplication-tables-build-') as directory:
        temporary = Path(directory)
        command = [compiler, '--keep-logs', '--keep-intermediates',
                   '--outdir', str(temporary)]
        if args.offline:
            command.append('--only-cached')
        command.append(str(ROOT / 'main.tex'))
        result = subprocess.run(command, cwd=ROOT, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if result.returncode:
            print(result.stdout, file=sys.stderr)
            return result.returncode
        pdf = temporary / 'main.pdf'
        bibliography = temporary / 'main.bbl'
        if not pdf.is_file() or not bibliography.is_file():
            print('The compiler did not produce the PDF and bibliography.', file=sys.stderr)
            return 1
        log = (temporary / 'main.log').read_text(errors='replace')
        if 'There were undefined references' in log or 'multiply defined' in log:
            print(log, file=sys.stderr)
            return 1
        shutil.copyfile(pdf, ROOT / PDF_NAME)
        for line in result.stdout.splitlines():
            if 'warning:' in line.lower():
                print(line, file=sys.stderr)
    print(f'Built {PDF_NAME}.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
