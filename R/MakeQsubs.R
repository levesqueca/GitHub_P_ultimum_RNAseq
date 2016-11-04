# Generating a function to make qsub jobs and a bash script to submit them:
# This makes a function that generates qsub commands and a bash to submit them in parallel on 
# processors on the biocluster using 3 inputs (will be used often throughout the script whenever a
# Linux command needs to be done on several data files)

MakeQsubs <- function(cmd, prefix, suffix = ".sub", node =1) {
  dir.create(paste(sharedPathAn, prefix, sep = ""), 
             showWarnings = TRUE, 
             recursive = FALSE)
  outPath <- paste(sharedPathAn, prefix, "/", sep="")
  for(k in 1:length(cmd)) {
    cat(paste("#!/bin/bash \n",
              "#$ -S /bin/bash
              # Ensure .e and .o and other output files go to working directory
              #$ -cwd
              # Request one slot in the smp environment
              #$ -pe smp ", node, "\n",
              "# Actual linux command for qsub \n",
              cmd[k],
              sep=""),
        file=paste(outPath, prefix, k, suffix, sep="")
    )
  }
  # make a bash script to run all qsub
  cat(paste("#!/bin/bash
            #$ -S /bin/bash
            argc=$#
            requiredArgc=0
            if [ $argc -ne $requiredArgc ]; then
            echo './test_mkdir.sh'
            exit 1
            fi
            
            prefixInFiles=", prefix, "\n",
            "suffixInFiles=", suffix, "\n",
            "for (( i = 1; i <= ", length(cmd), " ; i++ )); do 
            # keep track of what is going on...
            echo 'Treating file'  $prefixInFiles$i$suffixInFiles
            # define a script name that will be submited to the queue
            qsubFile=$prefixInFiles$i$suffixInFiles
            # make the script executable
            chmod a+x $qsubFile
            # submit the script to the queue
            qsub -cwd $qsubFile
            done", sep=""), 
      file=paste(outPath, prefix, ".sh", sep=""))
  cat(c("\n"," *** SUBMIT FOLLOWING TWO COMMANDS FROM HEADNODE ***","\n",
        "    1- Ensure working directory given for outputs:", "\n",
        paste("cd", outPath), "\n",
        "    2- Run the bash from within this working directory:", "\n",
        paste("bash ", prefix, ".sh", sep="")))
}
# End of qsub and bash making function