inpDir = params.inpDir
seqDir = params.seqDir
outDir = file(params.outDir)
sample = params.sample
def input = []
sample.eachWithIndex { it,i ->
	input[i] = file("${inpDir}/Input-$it*")
}
def fastq = []
sample.eachWithIndex { it,i ->
	fastq[i] = file("${seqDir}/$it*")
}
allFastQ = input + fastq
allFastQ.each { println it }
refDir = file(params.reference)
qualityScore = params.qualityScore
adapterSeq = file(params.adapterSeq)
genomebed = params.genomebed
outDir.mkdirs()
blacklist = file(params.blacklist)
println "Running pipeline for ${sample}"

process runFastQC{
	beforeScript { 'module load FastQC' }	
	publishDir "$outDir", mode: 'copy'

	input:
	each file('*') from allFastQ

	output:
	path "fastqc/" into fastqcpath
	file("fastqc/*.zip") into fastQC_files

	"""
	mkdir -p fastqc/
	fastqc --outdir fastqc/ -t1 *.fastq.gz
	"""
}

sampleNames = []
allFastQ.eachWithIndex { it,i ->
        sampleNames[i] = it.getAt(0).name.replaceAll(/R[1,2].fastq.gz/,"")
}

process align{
	tag "Aligning samples to GRCh38"
	cpus { 16 }
	memory { 32.GB }
	beforeScript {
		'module load bowtie2'
		'module load samtools'
	}
	publishDir "$outDir", mode: 'copy'

	input: 
		val sampleName from Channel.fromList(sampleNames)
		file('*') from Channel.fromList(allFastQ)
		path refDir
		
	output: 
		path "bam/" into bampath
		file("bam/${sampleName}raw.bam") into alignments

	"""
	mkdir -p bam
	bowtie2 -p 16 -q --very-sensitive-local -I 10 -X 700 -x $refDir/genome -1 ${sampleName}*R1* -2 ${sampleName}*R2* -S ${sampleName}.sam
	samtools view -bS ${sampleName}.sam > bam/${sampleName}raw.bam
	rm -rf ${sampleName}.sam
	"""
}

process dedup{
	cpus { 16 } 
	beforeScript { 
		'module load samtools'
		'module load bedtools'
	}
	memory { 64.GB }
	publishDir "$outDir/bam/", mode: 'copy'

	input:
		val sampleName from Channel.fromList(sampleNames)
		file("*") from alignments.collect()
		file blacklist

	output:
		file("${sampleName}raw.bam") 
		file("${sampleName}dedup.bam*") into dedup_aln
		file("${sampleName}dupmarked.bam") into dupmarked
		file("${sampleName}metrics.txt") into dupmark_metrics
		file("${sampleName}seq_metrics.csv")

	"""
	samtools view -bh -f 3 -F 4 -F 8 -q 20 ${sampleName}raw.bam | bedtools intersect -abam - -b $blacklist -v > ${sampleName}filtered.bam
	picard SortSam INPUT=${sampleName}filtered.bam OUTPUT=${sampleName}sorted.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
	picard MarkDuplicates INPUT=${sampleName}sorted.bam OUTPUT=${sampleName}dupmarked.bam VALIDATION_STRINGENCY=SILENT METRICS_FILE=${sampleName}metrics.txt
	samtools view -bh -F 1024 ${sampleName}dupmarked.bam > ${sampleName}dedup.bam
	samtools index ${sampleName}dedup.bam
	echo \$(samtools view -f 1 -c ${sampleName}raw.bam),\$(samtools view -f 3 -c ${sampleName}raw.bam),\$(samtools view -f 4 -c ${sampleName}raw.bam),,\$(samtools view -f 3 -F 4 -F 8 -F 512 -c ${sampleName}filtered.bam),\$(samtools view -f 1024 -c ${sampleName}dupmarked.bam) > ${sampleName}seq_metrics.csv
	"""			
}

dedup_aln.into { dedup_macs2; dedup_normBW }

process narrow_macs2{
	
	cpus { 16 }
	memory { 32.GB }
	publishDir "$outDir", mode: 'copy'

	input:
	val sampleName from Channel.fromList(sample)
	file '*' from dedup_macs2.collect()

	output:
	path "macs2/narrow/" into macs2path
	file("macs2/narrow/*")

	"""
	mkdir -p macs2/narrow
	source activate macs2
	macs2 callpeak -t ${sampleName}_dedup.bam -c Input-${sampleName}_dedup.bam -f BAMPE -n $sampleName --outdir macs2/narrow/ -g 1.87e9 --q 0.01 --broad --nomodel --keep-dup all
	"""
}

process normBW{

	cpus { 16 } 
	memory { 32.GB }
	
	publishDir "$outDir", mode: 'copy'
	

	input: 
	val sampleName from Channel.fromList(sampleNames)
	file('*') from dedup_normBW.collect()
	

	output:
	path "BWfiles" into bwfilespath
	file("BWfiles/${sampleName}rpkm.norm.bw")

	"""
	mkdir -p BWfiles
	source activate deeptools
	bamCoverage --bam ${sampleName}dedup.bam -o BWfiles/${sampleName}rpkm.norm.bw --binSize 10 --normalizeUsing RPKM --ignoreForNormalization chrX --skipNAs --extendReads 150 --numberOfProcessors 16 
	"""	

}

