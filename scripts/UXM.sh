#!/bin/bash

output=$1
input$2
atlas=$3
export PATH=${{PATH}}:/data/NHLBI_IDSS/references/UXM

module load samtools bedtools bamtools

uxm deconv ${input} -o ${output} --atlas ${atlas} --ignore Colon-Fibro Dermal-Fibro Gallbladder Bone-Osteob
