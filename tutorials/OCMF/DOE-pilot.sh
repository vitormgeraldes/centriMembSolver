#!/bin/bash
# create_DOE8.sh
# Run this script in the directory where OCMF_template/ is located

if [ -z "$BASH_VERSION" ]; then
    echo "Please run this script with bash, not sh."
    exit 1
fi

TEMPLATE="OCMF_template"
ROOT="OCMF_DOE8"

mkdir -p "$ROOT" || exit 1

make_case() {
    solute="$1"
    Pbar="$2"      # pressure [bar]
    Pref="$3"      # pressure [Pa]
    RPM="$4"
    Jlabel="$5"    # flux in LMH (for folder name)
    Jval="$6"      # flux in m/s (for Jv)
    C="$7"         # concentration [kg/m3]
    Cname="$8"     # for folder name

    caseName="${solute}_P${Pbar}bar_J${Jlabel}LMH_C${Cname}kgm3"
    echo "Creating: $caseName"

    cp -r "$TEMPLATE" "$ROOT/$caseName" || exit 1
    simFile="$ROOT/$caseName/constant/simulationConditions"

    # Update Jv, CA0, Pref, RPM
    sed -i \
        -e "s/^[[:space:]]*Jv[[:space:]].*/Jv          ${Jval};/" \
        -e "s/^[[:space:]]*CA0[[:space:]].*/CA0         ${C};/" \
        -e "s/^[[:space:]]*Pref[[:space:]].*/Pref        ${Pref};/" \
        -e "s/^[[:space:]]*RPM[[:space:]].*/RPM         ${RPM};/" \
        "$simFile"

    # Add header comment
    sed -i "1i // Case: $caseName" "$simFile"
}

# Flux levels (LMH → m/s)
J20=5.56e-06     # 20 LMH
J60=1.67e-05     # 60 LMH

# Pressures (your geometry)
# 10 bar -> 7434 rpm ; 30 bar -> 12876 rpm

########################
# NaCl – 4 cases
# C = 1 and 5 kg/m3
########################

# 1) NaCl, C=1, P=10 bar, J=20 LMH
make_case NaCl 10 1.0e6  7434 20 "$J20" 1.0 1
# 2) NaCl, C=5, P=10 bar, J=20 LMH
make_case NaCl 10 1.0e6  7434 20 "$J20" 5.0 5
# 3) NaCl, C=1, P=30 bar, J=60 LMH
make_case NaCl 30 3.0e6 12876 60 "$J60" 1.0 1
# 4) NaCl, C=5, P=30 bar, J=60 LMH
make_case NaCl 30 3.0e6 12876 60 "$J60" 5.0 5

########################
# Lactose – 4 cases
# C = 1 and 5 kg/m3
########################

# 5) Lactose, C=1, P=10 bar, J=20 LMH
make_case Lactose 10 1.0e6  7434 20 "$J20" 1.0 1
# 6) Lactose, C=5, P=10 bar, J=20 LMH
make_case Lactose 10 1.0e6  7434 20 "$J20" 5.0 5
# 7) Lactose, C=1, P=30 bar, J=60 LMH
make_case Lactose 30 3.0e6 12876 60 "$J60" 1.0 1
# 8) Lactose, C=5, P=30 bar, J=60 LMH
make_case Lactose 30 3.0e6 12876 60 "$J60" 5.0 5

echo "Done. DOE 8 cases created in: $ROOT/"
