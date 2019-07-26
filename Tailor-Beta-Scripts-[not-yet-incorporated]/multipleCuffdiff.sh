#!/bin/bash

#--------------------------------------------------------------------------------------------------------------
# multipleCuffdiff -Runs a differential gene expression analysis leaving out one sample each time in order to find significant changes in transcript expression, splicing, and promoter use that are consistently differentially expressed independent of biases from each individual sample.
#
#   - Run Cuffdiff by using the merged transcriptome assembly along with the BAM files from TopHat for each replicate
#   - Can send you emails when each job starts to execute and when each is finished so that you know when to submit the jobs for the next step
#--------------------------------------------------------------------------------------------------------------

#----SETTINGS---------------
source "${TAILOR_CONFIG}"
#---------------------------

#----COMMAND LINE ARGUMENTS-----------------------------------------------------------
readarray -t experiment1Group < $1   # read in all replicates, -t strips newlines
experiment1Name=$2

readarray -t experiment2Group < $3   # read in all replicates, -t strips newlines
experiment2Name=$4
#-------------------------------------------------------------------------------------
#----JOB SUBMISSION PARAMETERS---------------------------------------------------------------------
PROCESSORS=3
MEMORY="60536"      # PER PROCESSOR!! - 2048=2G, 4096=4G, 8192=8G, 16384=16G, 32768=32G, 65536=64G
DURATION="192:00"   # HH:MM - 72:00=3 days, 96:00=4 days, 192:00=8 days, 384:00=16 days, 768:00=32 days
QUEUE="long"        # short = max 4 hours;   long = max 30 days
#--------------------------------------------------------------------------------------------------

#----PATHS-----------------------------------------------------------------------------------------
INPUT_GTF="${CUFFDIFF_INPUT_GTF}"
INPUT_ALIGNMENTS="${LEAVEOUT_INPUT_ALIGNMENTS}"
OUTPUT="${LEAVEOUT_OUTPUT}"
SCRIPTS=${OUTPUT}/${JOBS_SCRIPTS_DIR}
JOBS=${OUTPUT}/${JOBS_OUT_DIR}
#--------------------------------------------------------------------------------------------------
#----OUTPUT------------------                                                                                                                    
if [ ! -d ${OUTPUT} ]; then
    mkdir ${OUTPUT}
fi
if [ ! -d ${SCRIPTS} ]; then
    mkdir ${SCRIPTS}
fi
if [ ! -d ${JOBS} ]; then
    mkdir ${JOBS}
fi
#----------------------------                    

#----GLOB INPUT--------------------------------
 findResults=${OUTPUT}/input.list.txt
 find ${INPUT_ALIGNMENTS} -type d -maxdepth 1 -name "${LEAVEOUT_FIND_INPUT_PATTERN}" > $findResults  
# search for only directories and no subdirs
#----------------------------------------------

COMMAND=cuffdiff
FOLD_CHANGE=$(eval echo $CUFFDIFF_FOLDCHANGE_OUTPUT)

#-------------------------------------------------------------------
# Create a list of all the BAM files for the 2 experiments
#-------------------------------------------------------------------

echo -e "\nGroup 1 size is ${#experiment1Group[@]}"
echo -e "\nGroup 1 is (${experiment1Group[@]})"
echo -e "\nGroup 2 size is ${#experiment2Group[@]}"
echo -e "\nGroup 2 is (${experiment2Group[@]})"

experiment1Files=""
experiment2Files=""

# ${!array[*]} gives indicies
# ${#array[@]} gives the length

while IFS='' read -r line || [[ -n "$line" ]]; do   # each line is a directory that contains a accepted_hit.bam file                             

    if [[ $line =~ $LEAVEOUT_REGEX ]]; then

        i=$(basename ${BASH_REMATCH[1]})
 
	LeaveOutFile=$i
    fi
    for j in ${!experiment1Group[*]}
    do
    # append a comma if necessary
	if [ $j -gt 0 ] && [ $j != $LeaveOutFile ]
	then
	    experiment1Files+=","
	fi
	experiment1Files+="${INPUT_ALIGNMENTS}/${experiment1Group[$j]}_out/${LEAVEOUT_ALIGNMENT_FILE}"
    done


    for j in ${!experiment2Group[*]}
    do
    # append a comma if necessary
	if [ $j -gt 0 ] && [ $j != $LeaveOutFile ]
	then
	    experiment2Files+=","
	fi
	experiment2Files+="${INPUT_ALIGNMENTS}/${experiment2Group[$j]}_out/${LEAVEOUT_ALIGNMENT_FILE}"
    done

#--------------------------------------------------------------------
#--------------------------------------------------------------------


# COMMAND_LINE="${COMMAND} -p ${PROCESSORS} -N -u -L ${experiment1Name},${experiment2Name} -library-norm-method classic-fpkm -b ${MULTIFASTA} -o ${OUTPUT}/${FOLD_CHANGE} ${INPUT_GTF}/merged.gtf $experiment1Files $experiment2Files"

# COMMAND_LINE="${COMMAND} -p ${PROCESSORS} -N -u -L ${experiment1Name},${experiment2Name} -dispersion-method per-condition -b ${MULTIFASTA} -o ${OUTPUT}/${FOLD_CHANGE} ${INPUT_GTF}/merged.gtf $experiment1Files $experiment2Files"

    EXTRA_PARAMETERS=$(eval echo "$EXTRA_LEAVEOUT_PARAMETERS")
    echo -e "\n${EXTRA_PARAMETERS}"

    COMMAND_LINE="${COMMAND} -p ${PROCESSORS} ${EXTRA_PARAMETERS} -o ${OUTPUT}/Leave_${i}_out/${FOLD_CHANGE} ${INPUT_GTF} $experiment1Files $experiment2Files"

    scriptString="mktemp -p ${SCRIPTS} -t ${COMMAND}.${i}.${FOLD_CHANGE}.XXXXXXXXXXX"
    echo -e "\n${scriptString}"
    tempScript=`${scriptString}`
    echo -e "\n${tempScript}"
    chmod=`chmod 777 ${tempScript}`
    chmodString="chmod 777 ${tempScript}"
    echo -e `${chmodString}`

    echo -e "source loadModules.sh\n\n" > ${tempScript}
    echo "$COMMAND_LINE" >> ${tempScript}


    if [ $SCHEDULER == "sge" ]; then
	SUBMIT_COMMAND="qsub -q $QUEUE -cwd -S /bin/bash -N ${FOLD_CHANGE} -pe smp ${PROCESSORS} -l h_rt=${DURATION},s_rt=${DURATION},vf=${MEMORY} -m eas -M ${USER_EMAIL} ${tempScript}"
    else
	SUBMIT_COMMAND="bsub -q $QUEUE -J ${FOLD_CHANGE}.${COMMAND} -n ${PROCESSORS} -R model==Intel_EM64T -R span[hosts=1] -R rusage[mem=${MEMORY}] -W ${DURATION} -u ${LSB_MAILTO} -B -o ${JOBS}/${COMMAND}.${i}.${FOLD_CHANGE}.%J.out -e ${JOBS}/${COMMAND}.${i}.${FOLD_CHANGE}.%J.error bash ${tempScript}"
    fi

    date=`date`
    echo -e "\n# $date\n"      >> ${OUTPUT}/${COMMAND}.jobs.log
    echo -e "\n# Job Script\n" >> ${OUTPUT}/${COMMAND}.jobs.log
    cat ${tempScript}          >> ${OUTPUT}/${COMMAND}.jobs.log
    echo -e "\n# Job Submission\n${SUBMIT_COMMAND}\n" >> ${OUTPUT}/${COMMAND}.jobs.log
    echo -e "\n#-------------------------------------------------------------------------------------------------------" >> ${OUTPUT}/${COMMAND}.jobs.log

    echo `${SUBMIT_COMMAND}`

done < "$findResults"  # Sends the $findResults file as input to the while-loop
