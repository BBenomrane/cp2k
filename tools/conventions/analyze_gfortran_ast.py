#!/usr/bin/env python3

# author: Ole Schuett

import argparse
import re
from os import path


BANNED_STM = ("GOTO", "FORALL", "OPEN", "CLOSE", "STOP")
BANNED_CALL = ("cp_fm_gemm", "m_abort", "mp_abort")
USE_EXCEPTIONS = ("omp_lib", "omp_lib_kinds", "lapack")

# precompile regex
re_symbol = re.compile(r"^\s*symtree.* symbol: '([^']+)'.*$")
re_use = re.compile(r" USE-ASSOC\(([^)]+)\)")
re_conv = re.compile(
    r"__(real_[48]_r8|real_4_c8|cmplx1_4_r8_r8)\[\["
)  # ignores integers


# ======================================================================================
def process_log_file(fhandle):
    public_symbols = set()
    used_symbols = set()

    def lprint(*args, **kwargs):
        return print("{}:".format(path.basename(fhandle.name)[:-4]), *args, **kwargs)

    module_name = None

    cur_sym = cur_proc = stat_var = stat_stm = None
    skip_until_DT_END = False

    for line in fhandle:
        line = line.strip()
        tokens = line.split()

        if skip_until_DT_END:
            # skip TRANSFERs which are part of READ/WRITE statement
            assert tokens[0] in ("DO", "TRANSFER", "END", "DT_END")
            if tokens[0] == "DT_END":
                skip_until_DT_END = False

        elif stat_var:
            if stat_var in line:  # Ok, something was done with stat_var
                stat_var = stat_stm = None  # reset
            elif line == "ENDIF":
                pass  # skip, check may happen in outer scope
            elif stat_stm == "ALLOCATE" and tokens[0] == "ASSIGN":
                pass  # skip lines, it's part of the ALLOCATE statement
            else:
                lprint(f'Found {stat_stm} with unchecked STAT in "{cur_proc}"')
                stat_var = stat_stm = None  # reset

        elif line.startswith("procedure name ="):
            cur_proc = line.split("=")[1].strip()
            if not module_name:
                module_name = cur_proc

        elif line.startswith("symtree: ") or len(line) == 0:
            cur_sym = None
            if len(line) == 0:
                continue
            cur_sym = re_symbol.match(line).group(1)

        elif line.startswith("attributes:"):
            is_imported = "USE-ASSOC" in line
            is_param = "PARAMETER" in line
            is_func = "FUNCTION" in line
            is_impl_save = "IMPLICIT-SAVE" in line
            is_impl_type = "IMPLICIT-TYPE" in line
            is_module_name = cur_proc == module_name

            if is_imported:
                mod = re_use.search(line).group(1)
                used_symbols.add(mod + "::" + cur_sym)
                if "MODULE  USE-ASSOC" in line and mod.lower() not in USE_EXCEPTIONS:
                    lprint(f'Module "{mod}" USEd without ONLY clause or not PRIVATE')

            # if(("SAVE" in line) and ("PARAMETER" not in line) and ("PUBLIC" in line)):
            #    print(loc+': Symbol "'+cur_sym+'" in procedure "'+cur_proc+'" is PUBLIC-SAVE')

            if is_impl_save and not is_param and not is_imported and not is_module_name:
                lprint(f'Symbol "{cur_sym}" in procedure "{cur_proc}" is IMPLICIT-SAVE')

            if is_impl_type and not is_imported and not is_func:
                lprint(f'Symbol "{cur_sym}" in procedure "{cur_proc}" is IMPLICIT-TYPE')

            if "THREADPRIVATE" in line:
                lprint(f'Symbol "{cur_sym}" in procedure "{cur_proc}" is THREADPRIVATE')

            if "PUBLIC" in line:
                public_symbols.add(module_name + "::" + cur_sym)

        elif line.startswith("!$OMP PARALLEL"):
            if "DEFAULT(NONE)" not in line:
                lprint(f'OMP PARALLEL without DEFAULT(NONE) found in "{cur_proc}"')

        elif line.startswith("CALL"):
            if tokens[1].lower() in BANNED_CALL:
                lprint(f'Found CALL {tokens[1]} in procedure "{cur_proc}"')
            elif tokens[1].lower().startswith("_gfortran_arandom_"):
                lprint(f'Found CALL RANDOM_NUMBER in procedure "{cur_proc}"')
            elif tokens[1].lower().startswith("_gfortran_random_seed_"):
                lprint(f'Found CALL RANDOM_SEED in procedure "{cur_proc}"')

        elif tokens and tokens[0] in BANNED_STM:
            lprint(f'Found {tokens[0]} statement in procedure "{cur_proc}"')

        elif line.startswith("WRITE"):
            unit = tokens[1].split("=")[1]
            if unit.isdigit():
                lprint(f'Found WRITE statement with hardcoded unit in "{cur_proc}"')

        elif line.startswith("DEALLOCATE") and "STAT=" in line:
            if ":ignore __final_" not in line:  # skip over auto-generated destructors
                lprint(f'Found DEALLOCATE with STAT argument in "{cur_proc}"')

        elif "STAT=" in line:  # catches also IOSTAT
            stat_var = line.split("STAT=", 1)[1].split()[0]
            stat_stm = line.split()[0]
            skip_until_DT_END = stat_stm in ("READ", "WRITE")

        elif "_gfortran_float" in line:
            lprint(f'Found FLOAT in "{cur_proc}"')

        elif re_conv.search(line):
            for m in re_conv.finditer(line):
                args = parse_args(line[m.end() :])
                if not re.match(r"\((kind = )?[48]\)", args[-1]):
                    lprint(
                        f'Found lossy conversion {m.group(1)} without KIND argument in "{cur_proc}"'
                    )

    # check for run-away DT_END search
    assert skip_until_DT_END is False

    return (public_symbols, used_symbols)


# ======================================================================================
def parse_args(line):
    assert line[0] == "("
    parentheses = 1
    args = list()
    for i in range(1, len(line)):
        if line[i] == "(":
            if parentheses == 1:
                a = i  # beginning of argument
            parentheses += 1
        elif line[i] == ")":
            parentheses -= 1
            if parentheses == 1:  # end of argument
                args.append(line[a : i + 1])
            if parentheses == 0:
                return args

    raise Exception("Could not find matching parentheses")


# ======================================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Checks the given ASTs for violations of the coding conventions",
        epilog="""\
For generating the abstract syntax tree (ast) run gfortran
with "-fdump-fortran-original" and redirect output to file.
This can be achieved by putting
    FCLOGPIPE = >$(notdir $<).ast
in the cp2k arch-file.
""",
    )
    parser.add_argument(
        "files",
        metavar="<ast-file>",
        type=str,
        nargs="+",
        help="files containing dumps of the AST",
    )
    args = parser.parse_args()

    for fn in args.files:
        assert fn.endswith(".ast")

        with open(fn, encoding="utf8") as fhandle:
            process_log_file(fhandle)

# EOF
