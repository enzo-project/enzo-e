#!/usr/bin/python3

# this is a simple tool used to update generated parser-files
# -> the conventions are based off of a comment from the old scons build-system
# -> https://github.com/enzo-project/enzo-e/blob/dfe1fdc602c37a556f036c44322257a92fb614c2/src/Cello/SConscript#L135

import argparse
import os
import shutil
import subprocess
import tempfile

_SOURCE_FILES = ['parse.l', 'parse.y']

parser = argparse.ArgumentParser(
    prog = 'regenerate_parser',
    description = (
        'This program is used to regenerates the C source and header files '
        'that are used for parameter-parsing after parse.y or parse.l has '
        'been modified. NOTE: requires "flex" and "bison"'
    )
)

parser.add_argument(
    'source_file', choices = ['both'] + _SOURCE_FILES,
    help = ('specifies which file(s) to use as inputs')
)

parser.add_argument(
    '--skip-cleanup-if-err', action = 'store_true',
    help = 'Enable this flag to skip cleanup if something goes wrong'
)

def _copy_files(src_l, dest_dir, verbose = True):
    out = []
    for src in src_l:
        out.append(f'{dest_dir}/{os.path.basename(src)}')
        print(f"copying {src} into {dest_dir}")
        shutil.copyfile(src, out[-1])
    return out

def _call(command, cwd, verbose = True):
    if verbose:
        print(f'executing `{command}`\n-> from {cwd}')
    subprocess.run(command, check = True, cwd = cwd, shell = True)

if __name__ == '__main__':
    args = parser.parse_args()

    cello_source_dir = os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '../src/Cello'
    ))

    # get the list of files
    if args.source_file == 'both':
        source_file_names = _SOURCE_FILES
    else:
        source_file_names = [args.source_file]

    source_paths = []
    for name in source_file_names:
        path = os.path.join(cello_source_dir, name)
        if not os.path.isfile(path): # sanity check!
            raise RuntimeError(f"something is wrong! can't find {path}")
        source_paths.append(path)

    # create a temporary directory
    tmpdirname = tempfile.mkdtemp()
    any_errs = False
    try:
        print(f'created temporary directory: {tmpdirname}')
        
        # recreate the directory structure previously used for generating the
        # parser-files (this slightly affects structure of output files)
        
        relocate_dir = f'{tmpdirname}/build/Cello'
        print(f'creating {relocate_dir}')
        os.makedirs(relocate_dir)
        for source_path in source_paths:
            _copy_files([source_path], relocate_dir)[0]
            fname = os.path.basename(source_path)

            if fname == 'parse.y':
                # the file-prefix arg only exists to match convention
                arg_l = ['-d', '-t',
                         '--file-prefix=build/Cello/parse',
                         'build/Cello/parse.y']
                _call(f'bison {" ".join(arg_l)}', cwd = tmpdirname)
                _copy_files(
                    [f'{tmpdirname}/build/Cello/{rslt}'
                     for rslt in ['parse.tab.c','parse.tab.h']],
                    cello_source_dir)
            else:
                _call('flex build/Cello/parse.l', cwd = tmpdirname)
                _copy_files([f'{tmpdirname}/lex.yy.c'], cello_source_dir)

    except:
        any_errs = True
        raise
    finally:
        if (any_errs and args.skip_cleanup_if_err):
            print('skipping cleanup')
        else:
            print(f'cleanup: deleting {tmpdirname}')
            shutil.rmtree(tmpdirname)

