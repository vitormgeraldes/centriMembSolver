#!/bin/bash
# create_DOE.sh
# Generates the full Design of Experiments (DOE) for constant permeate flux
# Run this script in the directory where OCMF_template/ is located

TEMPLATE="OCMF_template"
ROOT="DOE_cases"

mkdir -p "$ROOT" || exit 1

###############################################################################
# Factors of the DOE
###############################################################################

# Operating pressures [bar], internal representation [Pa], and corresponding RPM
# (RPM values based on your geometry and pressure–centrifugal equivalent)
pressures=(10 30 60)
pressures_Pa=(1.0e6 3.0e6 6.0e6)
rpms=(7434 12876 18210)

# Constant permeate flux levels (target average flux)
# Flux labels are used for folder naming
fluxLabels=(10 30 60)                   # LMH
fluxValues=(2.78e-06 8.33e-06 1.67e-05) # m/s

###############################################################################
# Function for creating each case folder
###############################################################################
create_case () {
    solute="$1"
    Pbar="$2"
    Pref="$3"
    RPM="$4"
    Jlabel="$5" # LMH
    Jval="$6"   # m/s
    C="$7"      # concentration [kg/m3]
    Cname="$8"  # string for folder name

    caseName="${solute}_P${Pbar}bar_J${Jlabel}LMH_C${Cname}kgm3"
    echo "Creating: $caseName"

    cp -r "$TEMPLATE" "$ROOT/$caseName" || exit 1
    simFile="$ROOT/$caseName/constant/simulationConditions"

    # Update the simulationConditions file
    sed -i \
        -e "s/^[[:space:]]*Jv[[:space:]].*/Jv          ${Jval};/" \
        -e "s/^[[:space:]]*CA0[[:space:]].*/CA0         ${C};/" \
        -e "s/^[[:space:]]*Pref[[:space:]].*/Pref        ${Pref};/" \
        -e "s/^[[:space:]]*RPM[[:space:]].*/RPM         ${RPM};/" \
        "$simFile"

    # Add header comment for traceability
    sed -i "1i // Case: $caseName" "$simFile"
}

###############################################################################
# SOLUTES & CONCENTRATION LEVELS
###############################################################################

# NaCl: 0.1, 1, 5 kg/m3
NaCl_concs=(0.1 1.0 5.0)
NaCl_names=(0p1 1 5)

# MgSO4: 0.1, 3, 10 kg/m3
Mg_concs=(0.1 3.0 10.0)
Mg_names=(0p1 3 10)

# LiCl: 0.1, 1, 3 kg/m3
Li_concs=(0.1 1.0 3.0)
Li_names=(0p1 1 3)

# Lactose: 1, 10, 30 kg/m3
Lac_concs=(1.0 10.0 30.0)
Lac_names=(1 10 30)

###############################################################################
# LOOP — full DOE generation
###############################################################################

# Generic looping structure for each solute
run_solute () {
    solute="$1"
    concArray=("${!2}")
    nameArray=("${!3}")

    for i in "${!pressures[@]}"; do
        Pbar="${pressures[$i]}"
        Pref="${pressures_Pa[$i]}"
        RPM="${rpms[$i]}"

        for j in "${!fluxLabels[@]}"; do
            Jlabel="${fluxLabels[$j]}"
            Jval="${fluxValues[$j]}"

            for k in "${!concArray[@]}"; do
                C="${concArray[$k]}"
                Cname="${nameArray[$k]}"
                create_case "$solute" "$Pbar" "$Pref" "$RPM" "$Jlabel" "$Jval" "$C" "$Cname"
            done
        done
    done
}

# Run for each solute:
run_solute "NaCl" NaCl_concs[@] NaCl_names[@]
run_solute "MgSO4" Mg_concs[@] Mg_names[@]
run_solute "LiCl" Li_concs[@] Li_names[@]
run_solute "Lactose" Lac_concs[@] Lac_names[@]

echo "All cases successfully generated in: $ROOT/"
