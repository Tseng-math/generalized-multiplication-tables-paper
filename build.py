#!/usr/bin/env python3
"""Build the standalone manuscript using Tectonic or latexmk."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
PDF_NAME = 'dyadic-divisor-products-and-rectangular-multiplication-tables.pdf'


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', help='Tectonic or latexmk executable or path')
    parser.add_argument('--offline', action='store_true',
                        help='Use only cached Tectonic resources')
    parser.add_argument('--bundle', help='Local Tectonic resource bundle')
    parser.add_argument('--cache-dir', help='Tectonic resource cache directory')
    args = parser.parse_args()

    choice = args.compiler or shutil.which('tectonic') or shutil.which('latexmk')
    if not choice:
        parser.error('Install Tectonic or latexmk, or specify --compiler.')
    executable = shutil.which(choice)
    if not executable:
        candidate = Path(choice).expanduser().resolve()
        if candidate.is_file() and os.access(candidate, os.X_OK):
            executable = str(candidate)
    if not executable:
        parser.error(f'Compiler not found: {choice}')
    executable = str(Path(executable).resolve())
    is_latexmk = 'latexmk' in Path(executable).name.lower()
    if is_latexmk and (args.offline or args.bundle or args.cache_dir):
        parser.error('--offline, --bundle and --cache-dir are Tectonic options.')

    env = os.environ.copy()
    if args.cache_dir:
        env['TECTONIC_CACHE_DIR'] = str(Path(args.cache_dir).expanduser().resolve())
    bundle = str(Path(args.bundle).expanduser().resolve()) if args.bundle else None
    with tempfile.TemporaryDirectory(prefix='tseng-paper-') as directory:
        build = Path(directory)
        shutil.copyfile(ROOT / 'main.tex', build / 'main.tex')
        if is_latexmk:
            command = [executable, '-pdf', '-interaction=nonstopmode',
                       '-halt-on-error', 'main.tex']
        else:
            command = [executable, '--keep-logs', '--keep-intermediates',
                       '--outdir', str(build)]
            if bundle:
                command += ['--bundle', bundle]
            if args.offline:
                command += ['--only-cached']
            command.append(str(build / 'main.tex'))
        result = subprocess.run(command, cwd=build, env=env,
                                capture_output=True, text=True)
        transcript = result.stdout + result.stderr
        log = transcript
        if (build / 'main.log').is_file():
            # Earlier passes can report references which the final pass
            # resolves. Check the final TeX log, not the multipass transcript.
            log = (build / 'main.log').read_text(errors='replace')
        diagnostics = re.findall(
            r'^(?:!|(?:warning|error):|LaTeX Warning|Package \S+ Warning|'
            r'Overfull \\[hv]box|Underfull \\[hv]box|.*undefined references|'
            r'.*undefined citations|.*multiply defined).*$'
            , log, re.I | re.M)
        output = build / 'main.pdf'
        if result.returncode or diagnostics or not output.is_file():
            print('Build failed or produced unresolved diagnostics.', file=sys.stderr)
            print('\n'.join(diagnostics) or (transcript + '\n' + log)[-6000:],
                  file=sys.stderr)
            raise SystemExit(1)
        shutil.copyfile(output, ROOT / PDF_NAME)
    print(f'Built {PDF_NAME}')


if __name__ == '__main__':
    main()
