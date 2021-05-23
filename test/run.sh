#!/bin/bash

# ==================================================================================================
#      Setup
# ==================================================================================================

# Global settings
CORES=10

# Color the spectrum!
COLOR_RED="\033[31m"
COLOR_GREEN="\033[32m"
COLOR_END="\033[0m"

# Change to the parent directory of where this script is located, which is the main grenepipe
# directory, so that we can call this script from anywhere, and run everything from there.
cd `dirname ${0}`/..
BASEPATH=`pwd -P`

# Remove old output
# ./test/clean.sh

# Copy the reference genome from the exampl dir to here, so that we can run the prepare step
# without cluttering the example directory.
if [[ ! -d ./test/reference ]]; then
    mkdir -p ./test/reference/
    cp ./example/TAIR10_chr_all.fa.gz ./test/reference/TAIR10_chr_all.fa.gz
fi

# Copy the samples table, so that we can change the paths without changing the original,
# We need to use a different sed separator here that cannot occur in paths, to avoid conflict.
cat ./test/samples_template.tsv | sed "s?#BASEPATH#?${BASEPATH}?g" > ./test/samples.tsv

# For the config, we need a bit of a special setup, as we need to replace some paths in
# different ways. We do this a bit more complicated setup here, in order to be able to work
# with the original config file, so that we do not need to keep track of copies and changes.
make_config() {
    local TARGET=$1

    cp ./config.yaml ${TARGET}
    sed -i "s?/path/to/data/samples.tsv?${BASEPATH}/test/samples.tsv?g" ${TARGET}
    sed -i "s?/path/to/data/genome.fa?${BASEPATH}/test/reference/TAIR10_chr_all.fa?g" ${TARGET}
    # cat ./test/config_template.yaml | sed "s?#BASEPATH#?${BASEPATH}?g" > ./test/config.yaml

    # Need an extra replacement step for threads. Might change in the future - this is a bit
    # volatile. But works for now.
    sed -i "s/threads: 12/threads: 6/g" ${TARGET}
}

# Make a config that we just use for the prep step. Makes it simpler than re-using one
# of the actual snakemake run config files that we produce later on.
make_config ./test/config.yaml

# ==================================================================================================
#      Run and Monitor Snakemake
# ==================================================================================================

# Count how many snakemake things we ran.
PASSCOUNT=0
FAILCOUNT=0

run_snakemake() {
    local CASE=$1
    local DIRECTORY=$2
    local EXTRA=$3

    # User output
    echo "[========" `date "+%F %T"` "========]"
    printf "${COLOR_GREEN}Running ${CASE}${COLOR_END}\n"

    # Run snakemake in the background.
    # Importantly, specify the conda prefix, so that the tools do not have to be loaded each time.
    snakemake \
        --use-conda \
        --conda-prefix ${BASEPATH}/test/conda-envs \
        --cores ${CORES} \
        --directory ${DIRECTORY} \
        ${EXTRA} \
        >> ${DIRECTORY}/${CASE}.log 2>&1 &
    PROC_ID=$!

    # Use the process ID to keep looping here while it is running,
    # and print some nice user output about the progress. This might not catch all progress,
    # in cases where multiple steps finish while we sleep here, but that is okay. It's for our
    # test cases only anyway, and we can live with that.
    PROGRESS=""
    while kill -0 "$PROC_ID" >/dev/null 2>&1; do
        sleep 1
        CURRENT=`egrep "[0-9]* of [0-9]* steps \([0-9.]*%\) done" ${DIRECTORY}/${CASE}.log | tail -n 1`
        if [[ "$CURRENT" != "$PROGRESS" ]]; then
            echo "    ${CURRENT}"
            PROGRESS=${CURRENT}
        fi
    done

    # Final user output for the test case
    SUCCESS=`grep "[0-9]* of [0-9]* steps ([0-9.]*%) done" ${DIRECTORY}/${CASE}.log | grep "100%"`
    NOTHING=`grep "Nothing to be done." ${DIRECTORY}/${CASE}.log`
    if [[ ! -z "$NOTHING" ]] ; then
        PASSCOUNT=$((PASSCOUNT+1))
        printf "${COLOR_GREEN}    Nothing to be done${COLOR_END}\n"
    elif [[ ! -z "$SUCCESS" ]]; then
        PASSCOUNT=$((PASSCOUNT+1))
        printf "${COLOR_GREEN}    Finished successfully${COLOR_END}\n"
    else
        FAILCOUNT=$((FAILCOUNT+1))
        printf "${COLOR_RED}    Error occurred${COLOR_END}\n"
        return 1
    fi
}

# ==================================================================================================
#      Prep Step
# ==================================================================================================

# Run the prepare step
run_snakemake "prep" "./test/" "--snakefile ./rules/prep.smk"
RESULT=$?
if [[ ${RESULT} != 0 ]]; then
    exit 1
fi

# Manual call of the prep step. Replaced by the above wrapper call.
# snakemake \
#     --use-conda \
#     --conda-prefix ${BASEPATH}/test/conda-envs \
#     --cores ${CORES} \
#     --directory ./test/ \
#     --snakefile ./rules/prep.smk \
#     &> ./test/prep.log

# ==================================================================================================
#      Test Cases
# ==================================================================================================

# Either get all scripts that we have, or the use provided ones via wildcard inclusion
DICTS=`ls ./test/cases/*.txt`
[[ "${1}" ]] && DICTS=`ls ./test/cases/*${1}*.txt`

# Now run all test cases. We always use the base config. This way, when we update the main
# grenepipe config file for new functionality, we only have to do minimal replacement here as well.
for DICT in ${DICTS} ; do

    # We collect all lines to be replaced from a "dictionary" script, one for each test,
    # and produce a new config file for this test, with all the lines replaced.
    # We also check that all the dictionary lines are actually in the original config,
    # in order to avoid that we accidentally changed our base config without noticing, which would mean
    # that we do not run the test the way we want, so better check.
    TARGET=$(basename "${DICT}" .txt)
    # echo "Running test case ${TARGET}"

    # Skip tests that start with an underscore.
    if [[ $TARGET = _* ]] ; then
        continue
    fi

    # Allow to re-run the test without doing everything again
    # (in particular, snakemake steps that already succceeded).
    if [[ ! -d ./test/out-${TARGET} ]]; then

        # Do the replacement
        mkdir ./test/out-${TARGET}
        make_config ./test/out-${TARGET}/config.yaml
        # cat ./test/config_template.yaml | sed "s?#BASEPATH#?${BASEPATH}?g" > ./test/out-${TARGET}/config.yaml
        while IFS="" read -r entry || [ -n "$p" ]
        do
            FROM=`echo "${entry}" | cut -f 1`
            TO=`echo "${entry}" | cut -f 2`
            # echo "from -${FROM}- to -${TO}-"

            # Test that the FROM part exists, to make sure that we are doing the right thing
            # if we replace config lines later on.
            if ! grep -q "${FROM}" ./test/out-${TARGET}/config.yaml; then
                printf "${COLOR_RED}Test setting does not exist: ${FROM}${COLOR_END}\n"
                exit 1
            fi

            sed -i "s/${FROM}/${TO}/g" ./test/out-${TARGET}/config.yaml
        done < ${DICT}
    fi

    # Now run snakemake on the new config file.
    run_snakemake "${TARGET}" "./test/out-${TARGET}"

    # Manual call of the test case. Replaced by the above wrapper call.
    # snakemake \
    #     --use-conda \
    #     --conda-prefix ${BASEPATH}/test/conda-envs \
    #     --cores ${CORES} \
    #     --directory ./test/out-${TARGET} \
    #     &> ./test/out-${TARGET}/snakemake.log
done

# Final user output
echo "[========" `date "+%F %T"` "========]"
printf "${COLOR_GREEN}PASS ${PASSCOUNT}${COLOR_END}\n"
printf "${COLOR_GREEN}FAIL ${FAILCOUNT}${COLOR_END}\n"
