#!/usr/bin/env bash

imageDir="../images"

mkdir -p $imageDir

for i in Singularity.* ; do 
    name=${i#*\.}
    containerName=$(echo -e "${name}.simg")
    echo -e "Building $containerName ..."
    singularity build ${imageDir}/${containerName} ${i}
    echo -e "DONE"
done
