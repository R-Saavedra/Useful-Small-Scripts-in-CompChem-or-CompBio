#!/usr/bin/env bash
#
# manual_counterpoise_input_generator.sh
# ---------------------------------------------------------------------------
# Builds the manual (non-automatic) set of Gaussian counterpoise inputs
# needed when Counterpoise=N cannot be used (e.g. because SCRF is present).
#
# For N fragments it writes, for every fragment i:
#   <name_i>.gjf                       fragment i alone, no ghosts
#   <name_i>_ghost<other names>.gjf    fragment i real + every OTHER fragment
#                                       as ghost (Bq) atoms, dimer/complex basis
# plus one:
#   COMPLETE_<header basename>.gjf     the full system, all atoms real
#
# All jobs are single points (any opt/freq in the header is stripped) at the
# geometry given in --coordinates. `nosymm` is added if not already present.
#
# ---------------------------------------------------------------------------
# USAGE
#   manual_counterpoise_input_generator.sh \
#       --coordinates FILE.xyz \
#       --fragments N "spec1" "spec2" ... END \
#       --names "name1" "name2" ... END \
#       --header TEMPLATE.gjh \
#       [--charges "q1" "q2" ... END] \
#       [--multiplicities "m1" "m2" ... END] \
#       [--total-charge Q] \
#       [--total-multiplicity M] \
#       [--outdir DIR] \
#       [--allow-unassigned]
#
# --coordinates FILE.xyz
#     Standard XYZ file: atom count, comment line, then "Element x y z" per
#     atom. Atom IDs used below are 1-based line numbers in that atom block
#     (i.e. the order they appear in the file), NOT reset per element.
#
# --fragments N "spec" ... END
#     N = number of fragments, followed by exactly N fragment specs,
#     terminated by the literal token END. Each spec is either:
#       - a comma-separated list of atom IDs and/or inclusive ID:ID ranges,
#         e.g. "1:32,43,56:66" or "33:42,44:54"
#       - the literal keyword "others" (case-insensitive), meaning "every
#         atom not explicitly claimed by any other fragment". At most one
#         fragment may use "others"; it is resolved after all explicit
#         specs are parsed. This is handy for the common 2-fragment case:
#         you only need to type the atom list for fragment 1.
#     Every atom in the coordinates file must end up in exactly one
#     fragment (use --allow-unassigned to permit leftover atoms instead;
#     not needed if you used "others").
#
# --names "name" ... END
#     Exactly N fragment names, same order as --fragments. Used to build
#     output filenames, so keep them filesystem-safe (letters/digits/-/_).
#
# --header TEMPLATE.gjh
#     A Gaussian header: the %-lines (%mem, %nproc, %chk, ...), the route
#     line ("#..."), and OPTIONALLY, in one of the next few non-blank lines
#     (a title line in between is fine and is skipped over), a charge/
#     multiplicity specification using Gaussian's own Counterpoise molecule-
#     spec convention:
#         Qtot,Mtot,Q1,M1,Q2,M2,...,QN,MN
#     (comma and/or space separated) -- the same line you would use in a
#     real `Counterpoise=N` job: total system first, then each fragment in
#     fragment order. Example header for the 2-fragment case (+1 complex,
#     +1 fragment 1, 0 fragment 2, all singlets):
#         %chk=job.chk
#         %mem=32GB
#         %nproc=24
#         # B3LYP/6-311G** SCRF=(CPCM,solvent=diethylether) counterpoise
#
#         1,1,1,1,0,1
#     Functional, basis set, dispersion, SCRF, grid, SCF convergence etc.
#     are all inherited from the route line verbatim (job-type keywords and
#     `counterpoise` are stripped automatically, see below). If a %chk= line
#     is present, it is rewritten per output file so jobs don't clobber each
#     other's checkpoints.
#
# --charges / --multiplicities (optional, override the header)
#     N values, same order as --fragments, terminated by END.
# --total-charge / --total-multiplicity (optional, override the header)
#     Charge/multiplicity for the COMPLETE (all-real) system.
#
#     Resolution order for every one of these four quantities, independently:
#       1. explicit CLI flag, if given
#       2. the header's charge/mult line, if present and valid
#       3. default: charge 0 / multiplicity 1 per fragment; total charge =
#          sum of fragment charges; total multiplicity = 1
#     Ghost jobs always use the REAL fragment's own charge/multiplicity --
#     never the sum with the ghosted fragment(s), since ghost atoms
#     contribute basis functions but no electrons.
#
# --outdir DIR (optional, default: current directory)
# --allow-unassigned (optional): don't error out on atoms belonging to no
#     fragment; they are simply left out of every file.
#
# ---------------------------------------------------------------------------
# EXAMPLE (2-fragment BSSE job, charges/mults taken from the header, second
#          fragment is "everything else")
#
#   ./manual_counterpoise_input_generator.sh \
#       --coordinates complex.xyz \
#       --fragments 2 "1:32,43,54:66" "others" END \
#       --names "FAD-ARG" "STY" END \
#       --header job.gjh
#
# produces: FAD-ARG.gjf, FAD-ARG_ghostSTY.gjf,
#           STY.gjf,      STY_ghostFAD-ARG.gjf,
#           COMPLETE_job.gjf
# ---------------------------------------------------------------------------

set -euo pipefail

if ((BASH_VERSINFO[0] < 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 3))); then
  echo "ERROR: this script needs bash >= 4.3 (namerefs/mapfile). Current: $BASH_VERSION" >&2
  exit 1
fi

print_usage() {
  local prog
  prog="$(basename "$0")"
  cat <<USAGE_EOF
$prog - generate the manual counterpoise (BSSE) Gaussian inputs

Use this when Counterpoise=N cannot be used (e.g. with SCRF/implicit solvent).
For N fragments it writes one single-point .gjf per job, all at the geometry
of the input .xyz:

  NAME.gjf               fragment alone, no ghost atoms
  NAME_ghostOTHER.gjf    fragment real + all other fragments as ghost (Bq)
                         atoms (OTHER = other names joined with '+')
  COMPLETE_HEADER.gjf    the whole system, all atoms real

USAGE
  $prog --coordinates FILE.xyz \\
      --fragments N "SPEC" ... END \\
      --names "NAME" ... END \\
      --header TEMPLATE.gjh \\
      [--charges "Q" ... END] [--multiplicities "M" ... END] \\
      [--total-charge Q] [--total-multiplicity M] \\
      [--outdir DIR] [--allow-unassigned]
  $prog -h | --help

REQUIRED
  --coordinates FILE.xyz
        XYZ file (atom count, comment line, then "Element x y z" per atom).
        Atom IDs are 1-based positions in the atom block.
  --fragments N "SPEC" ... END
        N fragments, each given as a comma-separated list of atom IDs and/or
        inclusive ID:ID ranges, e.g. "1:32,43,54:66". The keyword "others"
        (at most once) means every atom not claimed by another fragment.
        The list is closed by the literal token END. Every atom must end up
        in exactly one fragment.
  --names "NAME" ... END
        N names, same order as --fragments, used in the output filenames.
  --header TEMPLATE.gjh
        Gaussian header: %-lines (%mem, %nproc, %chk, ...) and the route line.
        Functional, basis, dispersion, SCRF, grid, SCF options are inherited
        as written. opt/freq/counterpoise are stripped, nosymm is added if
        missing, and %chk is rewritten per output file.
        Optionally, within the next few non-blank lines after the route
        (a title line in between is fine), a charge/multiplicity line in
        Gaussian's Counterpoise convention:
            Qtot,Mtot,Q1,M1,...,QN,MN         e.g.  1,1,0,1,1,1

OPTIONAL
  --charges "Q" ... END            N fragment charges (override the header)
  --multiplicities "M" ... END     N fragment multiplicities (override header)
  --total-charge Q                 charge of the COMPLETE system
  --total-multiplicity M           multiplicity of the COMPLETE system
        Resolution order for each: command-line flag > header line > default
        (fragment charge 0, multiplicity 1; total charge = sum of fragments;
        total multiplicity 1).
  --outdir DIR                     output directory (default: .)
  --allow-unassigned               drop atoms in no fragment instead of
                                   stopping with an error
  -h, --help                       show this help and exit

NOTES
  Ghost jobs use only the real fragment's own charge/multiplicity: ghost
  atoms carry basis functions but no nuclei and no electrons. For a closed-
  shell singlet, (sum of Z of the real atoms - charge) must be even.

EXAMPLE
  $prog \\
      --coordinates complex.xyz \\
      --fragments 2 "1:32,43,54:66" "others" END \\
      --names "FAD-ARG" "STY" END \\
      --header job.gjh

  writes FAD-ARG.gjf, FAD-ARG_ghostSTY.gjf, STY.gjf, STY_ghostFAD-ARG.gjf
  and COMPLETE_job.gjf.
USAGE_EOF
}

# --help / -h wins over everything else, wherever it appears
for _arg in "$@"; do
  case "$_arg" in
    -h|--help) print_usage; exit 0 ;;
  esac
done
if [[ $# -eq 0 ]]; then
  echo "$(basename "$0"): no arguments given. Try '$(basename "$0") --help'." >&2
  exit 1
fi

die()  { echo "ERROR: $*" >&2; exit 1; }
note() { echo "NOTE: $*" >&2; }
warn() { echo "WARNING: $*" >&2; }

# ------------------------------------------------------------------ args --
COORD_FILE=""
HEADER_FILE=""
N_FRAG=""
declare -a FRAG_SPECS=() FRAG_NAMES=()
declare -a CLI_FRAG_CHARGES=() CLI_FRAG_MULTS=()
CLI_TOTAL_CHARGE=""
CLI_TOTAL_MULT=""
CLI_TOTAL_MULT_SET=0
OUTDIR="."
ALLOW_UNASSIGNED=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --coordinates)
      [[ $# -ge 2 ]] || die "--coordinates needs a value"
      COORD_FILE="$2"; shift 2 ;;
    --header)
      [[ $# -ge 2 ]] || die "--header needs a value"
      HEADER_FILE="$2"; shift 2 ;;
    --fragments)
      shift
      [[ $# -ge 1 ]] || die "--fragments needs an integer count"
      N_FRAG="$1"; shift
      [[ "$N_FRAG" =~ ^[0-9]+$ ]] || die "--fragments count must be a positive integer, got '$N_FRAG'"
      FRAG_SPECS=()
      while [[ $# -gt 0 && "$1" != "END" ]]; do FRAG_SPECS+=("$1"); shift; done
      [[ "${1:-}" == "END" ]] || die "--fragments list must be terminated by END"
      shift
      ;;
    --names)
      shift
      FRAG_NAMES=()
      while [[ $# -gt 0 && "$1" != "END" ]]; do FRAG_NAMES+=("$1"); shift; done
      [[ "${1:-}" == "END" ]] || die "--names list must be terminated by END"
      shift
      ;;
    --charges)
      shift
      CLI_FRAG_CHARGES=()
      while [[ $# -gt 0 && "$1" != "END" ]]; do CLI_FRAG_CHARGES+=("$1"); shift; done
      [[ "${1:-}" == "END" ]] || die "--charges list must be terminated by END"
      shift
      ;;
    --multiplicities)
      shift
      CLI_FRAG_MULTS=()
      while [[ $# -gt 0 && "$1" != "END" ]]; do CLI_FRAG_MULTS+=("$1"); shift; done
      [[ "${1:-}" == "END" ]] || die "--multiplicities list must be terminated by END"
      shift
      ;;
    --total-charge)
      [[ $# -ge 2 ]] || die "--total-charge needs a value"
      CLI_TOTAL_CHARGE="$2"; shift 2 ;;
    --total-multiplicity)
      [[ $# -ge 2 ]] || die "--total-multiplicity needs a value"
      CLI_TOTAL_MULT="$2"; CLI_TOTAL_MULT_SET=1; shift 2 ;;
    --outdir)
      [[ $# -ge 2 ]] || die "--outdir needs a value"
      OUTDIR="$2"; shift 2 ;;
    --allow-unassigned)
      ALLOW_UNASSIGNED=1; shift ;;
    -h|--help)
      print_usage; exit 0 ;;
    *)
      die "unknown argument '$1' (see --help)" ;;
  esac
done

[[ -n "$COORD_FILE" ]] || die "--coordinates is required"
[[ -f "$COORD_FILE" ]] || die "coordinates file '$COORD_FILE' not found"
[[ -n "$HEADER_FILE" ]] || die "--header is required"
[[ -f "$HEADER_FILE" ]] || die "header file '$HEADER_FILE' not found"
[[ -n "$N_FRAG" ]] || die "--fragments is required"
(( N_FRAG >= 1 )) || die "--fragments count must be >= 1"
[[ "${#FRAG_SPECS[@]}" -eq "$N_FRAG" ]] || die "expected $N_FRAG fragment specs, got ${#FRAG_SPECS[@]}"
[[ "${#FRAG_NAMES[@]}" -eq "$N_FRAG" ]] || die "expected $N_FRAG fragment names, got ${#FRAG_NAMES[@]}"

mkdir -p "$OUTDIR"

# --------------------------------------------------------------- read xyz --
mapfile -t XYZ_LINES < "$COORD_FILE"
[[ "${#XYZ_LINES[@]}" -ge 2 ]] || die "coordinates file too short to be a valid xyz"
NATOMS="$(tr -d '[:space:]' <<< "${XYZ_LINES[0]}")"
[[ "$NATOMS" =~ ^[0-9]+$ ]] || die "first line of '$COORD_FILE' must be an atom count"
(( ${#XYZ_LINES[@]} >= NATOMS + 2 )) || die "'$COORD_FILE' declares $NATOMS atoms but has too few lines"

declare -a ELEM X Y Z
for ((i=1; i<=NATOMS; i++)); do
  line="${XYZ_LINES[$((i+1))]}"
  read -r el x y z _ <<< "$line"
  [[ -n "$el" && -n "$x" && -n "$y" && -n "$z" ]] || die "could not parse atom $i: '$line'"
  ELEM[i]="$el"; X[i]="$x"; Y[i]="$y"; Z[i]="$z"
done

# ------------------------------------------------------- fragment mapping --
parse_spec_into() {
  # $1 = spec string, $2 = name of array variable to fill with atom IDs
  local spec="$1"
  local -n _out="$2"
  _out=()
  local tok lo hi k
  local IFS=','
  for tok in $spec; do
    if [[ "$tok" == *:* ]]; then
      lo="${tok%%:*}"; hi="${tok##*:}"
      [[ "$lo" =~ ^[0-9]+$ && "$hi" =~ ^[0-9]+$ ]] || die "bad range '$tok' in spec '$spec'"
      (( lo <= hi )) || die "range '$tok' has start > end"
      for ((k=lo; k<=hi; k++)); do _out+=("$k"); done
    else
      [[ "$tok" =~ ^[0-9]+$ ]] || die "bad atom id '$tok' in spec '$spec'"
      _out+=("$tok")
    fi
  done
}

declare -A FRAG_OF=()          # atom id -> fragment index (1-based)
declare -A FRAG_IDS=()         # fragment index -> space-separated sorted atom ids
OTHERS_FIDX=""

# pass 1: explicit specs (skip any "others" fragment for now)
for ((fidx=1; fidx<=N_FRAG; fidx++)); do
  spec="${FRAG_SPECS[$((fidx-1))]}"
  spec_lc="${spec,,}"
  if [[ "$spec_lc" == "others" ]]; then
    [[ -z "$OTHERS_FIDX" ]] || die "only one fragment may use 'others' (already used by fragment $OTHERS_FIDX, '${FRAG_NAMES[$((OTHERS_FIDX-1))]}')"
    OTHERS_FIDX=$fidx
    continue
  fi
  parse_spec_into "$spec" ids
  for id in "${ids[@]}"; do
    (( id >= 1 && id <= NATOMS )) || die "atom id $id (fragment $fidx, '${FRAG_NAMES[$((fidx-1))]}') out of range 1..$NATOMS"
    if [[ -n "${FRAG_OF[$id]:-}" ]]; then
      die "atom id $id assigned to both fragment ${FRAG_OF[$id]} and fragment $fidx"
    fi
    FRAG_OF[$id]=$fidx
  done
  mapfile -t sorted_ids < <(printf '%s\n' "${ids[@]}" | sort -n)
  FRAG_IDS[$fidx]="${sorted_ids[*]}"
done

# pass 2: resolve "others", if used, to every currently-unclaimed atom
if [[ -n "$OTHERS_FIDX" ]]; then
  others_ids=()
  for ((i=1; i<=NATOMS; i++)); do
    [[ -n "${FRAG_OF[$i]:-}" ]] || others_ids+=("$i")
  done
  [[ "${#others_ids[@]}" -gt 0 ]] || die "'others' fragment ('${FRAG_NAMES[$((OTHERS_FIDX-1))]}') would contain 0 atoms -- every atom is already claimed"
  for id in "${others_ids[@]}"; do FRAG_OF[$id]=$OTHERS_FIDX; done
  FRAG_IDS[$OTHERS_FIDX]="${others_ids[*]}"
  note "fragment '${FRAG_NAMES[$((OTHERS_FIDX-1))]}' (others) resolved to ${#others_ids[@]} atom(s): ${others_ids[*]}"
fi

unassigned=()
for ((i=1; i<=NATOMS; i++)); do
  [[ -n "${FRAG_OF[$i]:-}" ]] || unassigned+=("$i")
done
if [[ "${#unassigned[@]}" -gt 0 ]]; then
  if [[ "$ALLOW_UNASSIGNED" -eq 1 ]]; then
    warn "${#unassigned[@]} atom(s) belong to no fragment and will be dropped from every file: ${unassigned[*]}"
  else
    die "${#unassigned[@]} atom(s) belong to no fragment: ${unassigned[*]} (use --allow-unassigned to drop them intentionally, or a 'others' fragment to catch them)"
  fi
fi

# ------------------------------------------------------------ read header --
declare -a LINK0=()
ROUTE=""
declare -a POST_LINES=()
in_route=0
route_done=0
while IFS= read -r line || [[ -n "$line" ]]; do
  line="${line%$'\r'}"   # tolerate CRLF-terminated (Windows-edited) header files
  if [[ $route_done -eq 1 ]]; then
    POST_LINES+=("$line")
    continue
  fi
  if [[ "$line" == %* ]]; then
    LINK0+=("$line")
  elif [[ "$line" == \#* ]]; then
    ROUTE="$line"; in_route=1
  elif [[ $in_route -eq 1 ]]; then
    if [[ -z "${line// /}" ]]; then
      route_done=1
    else
      ROUTE="$ROUTE $line"
    fi
  fi
done < "$HEADER_FILE"
[[ -n "$ROUTE" ]] || die "no route line (starting with #) found in '$HEADER_FILE'"

# Candidate charge/mult lines: every non-blank line after the route (up to a
# handful of them), not just the very first one -- a title line is often
# present between the route and the actual charge/mult line.
declare -a HEADER_CM_CANDIDATES=()
for l in "${POST_LINES[@]}"; do
  [[ -n "${l// /}" ]] && HEADER_CM_CANDIDATES+=("$l")
  [[ "${#HEADER_CM_CANDIDATES[@]}" -lt 6 ]] || break
done

CLEAN_ROUTE="$ROUTE"
CLEAN_ROUTE=$(sed -E 's/[[:space:]]+\<counterpoise\>(=[0-9]+)?//Ig'  <<< "$CLEAN_ROUTE")
CLEAN_ROUTE=$(sed -E 's/[[:space:]]+\<opt\>(\([^)]*\))?(=[^[:space:]]+)?//Ig' <<< "$CLEAN_ROUTE")
CLEAN_ROUTE=$(sed -E 's/[[:space:]]+\<freq\>(\([^)]*\))?(=[^[:space:]]+)?//Ig' <<< "$CLEAN_ROUTE")
CLEAN_ROUTE=$(tr -s ' ' <<< "$CLEAN_ROUTE")
if [[ "$ROUTE" != "$CLEAN_ROUTE" ]]; then
  note "stripped opt/freq/counterpoise keyword(s) from route (these jobs are forced single point, manual BSSE)"
fi
if ! grep -qiw nosymm <<< "$CLEAN_ROUTE"; then
  CLEAN_ROUTE="$CLEAN_ROUTE nosymm"
  note "added nosymm to the route (recommended for manual BSSE / ghost-atom jobs)"
fi

header_base="$(basename "${HEADER_FILE%.*}")"
system_title="$(basename "${COORD_FILE%.*}") manual CP"

# ------------------------------------------- parse header charge/mult line --
HDR_OK=0
HDR_TOTAL_CHARGE=""
HDR_TOTAL_MULT=""
declare -a HDR_FRAG_CHARGES=() HDR_FRAG_MULTS=()

expected=$(( 2 * (N_FRAG + 1) ))
warned_malformed=0
for cand in "${HEADER_CM_CANDIDATES[@]}"; do
  [[ "$HDR_OK" -eq 0 ]] || break
  norm="${cand//,/ }"
  read -ra hdr_toks <<< "$norm"
  all_int=1
  for t in "${hdr_toks[@]}"; do
    [[ "$t" =~ ^[+-]?[0-9]+$ ]] || { all_int=0; break; }
  done
  if [[ "$all_int" -eq 1 && "${#hdr_toks[@]}" -eq "$expected" ]]; then
    HDR_TOTAL_CHARGE="${hdr_toks[0]}"
    HDR_TOTAL_MULT="${hdr_toks[1]}"
    for ((i=0; i<N_FRAG; i++)); do
      HDR_FRAG_CHARGES+=("${hdr_toks[$((2*(i+1)))]}")
      HDR_FRAG_MULTS+=("${hdr_toks[$((2*(i+1)+1))]}")
    done
    HDR_OK=1
    note "found charge/multiplicity line in '$HEADER_FILE': '$cand'"
  elif [[ "$all_int" -eq 1 && "$warned_malformed" -eq 0 ]]; then
    warn "found a charge/multiplicity-looking line in '$HEADER_FILE' ('$cand') with ${#hdr_toks[@]} value(s), expected $expected (2 per fragment + 2 for the total, for $N_FRAG fragments); ignoring it and checking for another"
    warned_malformed=1
  fi
done

# --------------------------------------------- resolve charges/mults/totals --
declare -a FRAG_CHARGES=() FRAG_MULTS=()

if [[ "${#CLI_FRAG_CHARGES[@]}" -gt 0 ]]; then
  [[ "${#CLI_FRAG_CHARGES[@]}" -eq "$N_FRAG" ]] || die "expected $N_FRAG charges via --charges, got ${#CLI_FRAG_CHARGES[@]}"
  FRAG_CHARGES=("${CLI_FRAG_CHARGES[@]}")
  note "fragment charges taken from --charges"
elif [[ "$HDR_OK" -eq 1 ]]; then
  FRAG_CHARGES=("${HDR_FRAG_CHARGES[@]}")
  note "fragment charges taken from the charge/multiplicity line in '$HEADER_FILE'"
else
  for ((i=0;i<N_FRAG;i++)); do FRAG_CHARGES+=("0"); done
  note "no fragment charges given (--charges or a charge/mult line in the header); defaulting every fragment to charge 0"
fi

if [[ "${#CLI_FRAG_MULTS[@]}" -gt 0 ]]; then
  [[ "${#CLI_FRAG_MULTS[@]}" -eq "$N_FRAG" ]] || die "expected $N_FRAG multiplicities via --multiplicities, got ${#CLI_FRAG_MULTS[@]}"
  FRAG_MULTS=("${CLI_FRAG_MULTS[@]}")
  note "fragment multiplicities taken from --multiplicities"
elif [[ "$HDR_OK" -eq 1 ]]; then
  FRAG_MULTS=("${HDR_FRAG_MULTS[@]}")
  note "fragment multiplicities taken from the charge/multiplicity line in '$HEADER_FILE'"
else
  for ((i=0;i<N_FRAG;i++)); do FRAG_MULTS+=("1"); done
  note "no fragment multiplicities given; defaulting every fragment to multiplicity 1"
fi

if [[ -n "$CLI_TOTAL_CHARGE" ]]; then
  TOTAL_CHARGE="$CLI_TOTAL_CHARGE"
  note "total charge taken from --total-charge"
elif [[ "$HDR_OK" -eq 1 ]]; then
  TOTAL_CHARGE="$HDR_TOTAL_CHARGE"
  note "total charge taken from the charge/multiplicity line in '$HEADER_FILE'"
else
  s=0; for q in "${FRAG_CHARGES[@]}"; do s=$((s+q)); done
  TOTAL_CHARGE=$s
  note "no total charge given; using sum of fragment charges = $TOTAL_CHARGE"
fi

if [[ "$CLI_TOTAL_MULT_SET" -eq 1 ]]; then
  TOTAL_MULT="$CLI_TOTAL_MULT"
  note "total multiplicity taken from --total-multiplicity"
elif [[ "$HDR_OK" -eq 1 ]]; then
  TOTAL_MULT="$HDR_TOTAL_MULT"
  note "total multiplicity taken from the charge/multiplicity line in '$HEADER_FILE'"
else
  TOTAL_MULT=1
  for m in "${FRAG_MULTS[@]}"; do
    if [[ "$m" != "1" ]]; then
      warn "fragment multiplicities are not all singlets; defaulting total multiplicity to 1 -- verify or pass --total-multiplicity"
      break
    fi
  done
fi

# ---------------------------------------------------------------- writers --
write_link0() {
  local outbase="$1" l
  for l in "${LINK0[@]}"; do
    if [[ "$l" =~ ^%[Cc][Hh][Kk]= ]]; then
      echo "%chk=${outbase}.chk"
    else
      echo "$l"
    fi
  done
}

write_gjf() {
  # $1 outfile  $2 title  $3 charge  $4 mult  $5 real_ids arrayname  $6 ghost_ids arrayname
  local outfile="$1" title="$2" charge="$3" mult="$4"
  local -n real_ref="$5"
  local -n ghost_ref="$6"
  local outbase
  outbase="$(basename "${outfile%.gjf}")"
  {
    write_link0 "$outbase"
    echo "$CLEAN_ROUTE"
    echo
    echo "$title"
    echo
    echo "$charge $mult"
    local id
    for id in "${real_ref[@]}"; do
      printf '%-8s%22s%22s%22s\n' "${ELEM[id]}" "${X[id]}" "${Y[id]}" "${Z[id]}"
    done
    for id in "${ghost_ref[@]}"; do
      printf '%-8s%22s%22s%22s\n' "${ELEM[id]}-Bq" "${X[id]}" "${Y[id]}" "${Z[id]}"
    done
    echo
  } > "$outfile"
}

printf '\n%-34s %6s %6s %8s %6s\n' "file" "real" "ghost" "charge" "mult"
printf '%.0s-' {1..70}; echo

empty_ids=()

for ((fidx=1; fidx<=N_FRAG; fidx++)); do
  name="${FRAG_NAMES[$((fidx-1))]}"
  charge="${FRAG_CHARGES[$((fidx-1))]}"
  mult="${FRAG_MULTS[$((fidx-1))]}"
  read -ra real_ids <<< "${FRAG_IDS[$fidx]}"

  ghost_ids=()
  other_names=()
  for ((j=1; j<=N_FRAG; j++)); do
    if [[ $j -ne $fidx ]]; then
      read -ra tmp <<< "${FRAG_IDS[$j]}"
      ghost_ids+=("${tmp[@]}")
      other_names+=("${FRAG_NAMES[$((j-1))]}")
    fi
  done
  mapfile -t ghost_ids < <(printf '%s\n' "${ghost_ids[@]}" | sort -n)

  out1="${OUTDIR}/${name}.gjf"
  write_gjf "$out1" "${system_title} - ${name} alone" "$charge" "$mult" real_ids empty_ids
  printf '%-34s %6d %6d %8s %6s\n' "$(basename "$out1")" "${#real_ids[@]}" 0 "$charge" "$mult"

  ghost_label=$(IFS=+; echo "${other_names[*]}")
  out2="${OUTDIR}/${name}_ghost${ghost_label}.gjf"
  write_gjf "$out2" "${system_title} - ${name} + ghost(${other_names[*]})" "$charge" "$mult" real_ids ghost_ids
  printf '%-34s %6d %6d %8s %6s\n' "$(basename "$out2")" "${#real_ids[@]}" "${#ghost_ids[@]}" "$charge" "$mult"
done

all_ids=()
for ((i=1; i<=NATOMS; i++)); do
  [[ -n "${FRAG_OF[$i]:-}" ]] && all_ids+=("$i")
done
complete_out="${OUTDIR}/COMPLETE_${header_base}.gjf"
write_gjf "$complete_out" "${system_title} - complete system" "$TOTAL_CHARGE" "$TOTAL_MULT" all_ids empty_ids
printf '%-34s %6d %6d %8s %6s\n' "$(basename "$complete_out")" "${#all_ids[@]}" 0 "$TOTAL_CHARGE" "$TOTAL_MULT"

echo
echo "route used: $CLEAN_ROUTE"
echo "output dir: $OUTDIR"
