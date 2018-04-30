#those functions take arguments from configuration file and return the tool path and an array with configurations info

def Set_bamsurgeon(conf,opts):
	ms_args=[]
	ms_args+=['--picardjar',conf["picard"]["path"]]
	if conf["bamsurgeon"]["aligner"]!="":
		ms_args+=['--aligner',conf["bamsurgeon"]["aligner"]]
	else:
		ms_args+=['--aligner','mem']
	if conf["bamsurgeon"]["snvfrac"]!="":
		ms_args+=['--snvfrac',conf["bamsurgeon"]["snvfrac"]]
	else:
		pass
	if conf["bamsurgeon"]["mutfrac"]!="":
		ms_args+=['--mutfrac',conf["bamsurgeon"]["mutfrac"]]
	else:
		pass
	if conf["bamsurgeon"]["numsnvs"]!="":
		ms_args+=['--numsnvs',conf["bamsurgeon"]["numsnvs"]]
	else:
		pass
	if conf["bamsurgeon"]["coverdiff"]!="":
		ms_args+=['--coverdiff',conf["bamsurgeon"]["coverdiff"]]
	else:
		pass
	if conf["bamsurgeon"]["haplosize"]!="":
		ms_args+=['--haplosize',conf["bamsurgeon"]["haplosize"]]
	else:
		pass
	if conf["bamsurgeon"]["procs"]!="":
		ms_args+=['--procs', conf["bamsurgeon"]["procs"]]
	else:
		pass
	if conf["bamsurgeon"]["mindepth"]!="":
		ms_args+=['--mindepth',conf["bamsurgeon"]["mindepth"]]
	else:
		pass
	if conf["bamsurgeon"]["maxdepth"]!="":
		ms_args+=['--maxdepth',conf["bamsurgeon"]["maxdepth"]]
	else:
		pass
	if conf["bamsurgeon"]["minmutreads"]!="":
		ms_args+=['--minmutreads',conf["bamsurgeon"]["minmutreads"]]
	else:
		pass
	if conf["bamsurgeon"]["avoidreads"]!="":
		ms_args+=['--avoidreads',conf["bamsurgeon"]["avoidreads"]]
	else:
		pass
	if conf["bamsurgeon"]["maxopen"]!="":
		ms_args+=['--maxopen',conf["bamsurgeon"]["maxopen"]]
	else:
		pass
	if conf["bamsurgeon"]["seed"]!="":
		ms_args+=['--seed',conf["bamsurgeon"]["seed"]]
	else:
		pass
	if conf["bamsurgeon"]["flags"]!="":
		ms_args+=conf["bamsurgeon"]["flags"].split(',')
	else:
		pass
	path_addsnv=conf["bamsurgeon"]["addsnv"]
	path_addindel=conf["bamsurgeon"]["addindel"]

	return path_addsnv,path_addindel,ms_args

# def Set_bamsurgeon_sv(conf,opts):
# 	ms_args=[]
# 	if conf["bamsurgeon_sv"]["aligner"]!="":
# 		ms_args+=['--aligner',conf["bamsurgeon_sv"]["aligner"]]
# 	else:
# 		ms_args+=['--aligner','mem']
# 	if conf["bamsurgeon_sv"]["svfrac"]!="":
# 		ms_args+=['--mutfrac',conf["bamsurgeon_sv"]["svfrac"]]
# 	else:
# 		pass
# 	if conf["bamsurgeon_sv"]["maxlibsize"]!="":
# 		ms_args+=['--maxlibsize',conf["bamsurgeon_sv"]["maxlibsize"]]
# 	else:
# 		pass
# 	if conf["bamsurgeon_sv"]["kmer"]!="":
# 		ms_args+=['--kmer',conf["bamsurgeon_sv"]["kmer"]]
# 	else:
# 		pass
# 	if conf["bamsurgeon_sv"]["maxctglen"]!="":
# 		ms_args+=['--maxctglen',conf["bamsurgeon_sv"]["maxctglen"]]
# 	else:
# 		pass
# 	if conf["bamsurgeon_sv"]["maxmut"]!="":
# 		ms_args+=['-n',conf["bamsurgeon_sv"]["maxmut"]]
# 	else:
# 		pass
# 	if conf["bamsurgeon_sv"]["procs"]!="":
# 		ms_args+=['--procs', conf["bamsurgeon_sv"]["procs"]]
# 	else:
# 		pass
# 	if conf["bamsurgeon_sv"]["ismean"]!="":
# 		ms_args+=['--ismean',conf["bamsurgeon_sv"]["ismean"]]
# 	else:
# 		pass
# 	if conf["bamsurgeon_sv"]["issd"]!="":
# 		ms_args+=['--issd',conf["bamsurgeon_sv"]["issd"]]
# 	else:
# 		pass
# 	if conf["bamsurgeon_sv"]["delay"]!="":
# 		ms_args+=['--delay',conf["bamsurgeon_sv"]["delay"]]
# 	else:
# 		pass
# 	if conf["bamsurgeon_sv"]["cnvfile"]!="":
# 		ms_args+=['--cnvfile',conf["bamsurgeon_sv"]["cnvfile"]]
# 	else:
# 		pass
# 	if conf["bamsurgeon_sv"]["seed"]!="":
# 		ms_args+=['--seed',conf["bamsurgeon_sv"]["seed"]]
# 	else:
# 		pass
# 	if conf["bamsurgeon_sv"]["flags"]!="":
# 		ms_args+=conf["bamsurgeon_sv"]["flags"].split(',')
# 	else:
# 		pass
# 	path_addsv=conf["bamsurgeon_sv"]["addsv"]

# 	return path_addsv,ms_args

def Set_art(conf,opts):
	art_args=[]
	if conf["art"]["qprof1"]!="":
		art_args+= ['--qprof1', conf["art"]["qprof1"]]
	else:
		pass
	if conf["art"]["qprof2"]!="":
		art_args+= ['--qprof2', conf["art"]["qprof2"]]
	else:
		pass
	if conf["art"]["rcount"]!="":
		art_args+= ['--rcount', conf["art"]["rcount"]]
	else:
		pass
	if conf["art"]["id"]!="":
		art_args+= ['--id', conf["art"]["id"]]
	else:
		pass
	if conf["art"]["insRate"]!="":
		art_args+= ['--insRate', conf["art"]["insRate"]]
	else:
		pass
	if conf["art"]["insRate2"]!="":
		art_args+= ['--insRate2', conf["art"]["insRate2"]]
	else:
		pass
	if conf["art"]["delRate"]!="":
		art_args+= ['--delRate', conf["art"]["delRate"]]
	else:
		pass
	if conf["art"]["delRate2"]!="":
		art_args+= ['--delRate2', conf["art"]["delRate2"]]
	else:
		pass
	if conf["art"]["maxIndel"]!="":
		art_args+= ['--maxIndel', conf["art"]["maxIndel"]]
	else:
		pass
	if conf["art"]["len"]!="":
		art_args+= ['--len', conf["art"]["len"]]
	else:
		art_args+= ['--len','150']
	if conf["art"]["mflen"]!="":
		art_args+= ['--mflen', conf["art"]["mflen"]]
	else:
		art_args+= ['--mflen','300']
	if conf["art"]["maskN"]!="":
		art_args+= ['--maskN', conf["art"]["maskN"]]
	else:
		pass
	if conf["art"]["minQ"]!="":
		art_args+= ['--minQ', conf["art"]["minQ"]]
	else:
		pass
	if conf["art"]["maxQ"]!="":
		art_args+= ['--maxQ', conf["art"]["maxQ"]]
	else:
		pass
	if conf["art"]["qShift"]!="":
		art_args+= ['--qShift', conf["art"]["qShift"]]
	else:
		pass
	if conf["art"]["qShift2"]!="":
		art_args+= ['--qShift2', conf["art"]["qShift2"]]
	else:
		pass
	if conf["art"]["rndSeed"]!="":
		art_args+= ['--rndSeed', conf["art"]["rndSeed"]]
	else:
		pass
	if conf["art"]["fcov"]!="":
		art_args+= ['--fcov', conf["art"]["fcov"]]
	else:
		art_args+= ['--fcov', '1000']
	if conf["art"]["sdev"]!="":
		art_args+= ['--sdev', conf["art"]["sdev"]]
	else:
		art_args+= ['--sdev', '10']
	if conf["art"]["seqSys"]!="":
		art_args+= ['--seqSys', conf["art"]["seqSys"]]
	else:
		art_args+= ['--seqSys', 'MSv3']
	art_args+=conf["art"]["flags"].split(',')
	path_art=conf["art"]["path"]
	if conf["art"]["prefix"]!="":
		art_args+= ['--out', opts.out_path+'/FASTQ/' +conf["art"]["prefix"]]
	else:
		art_args+= ['--out', opts.out_path+'/FASTQ/' + opts.sample]
	return path_art,art_args



def Set_bwa(conf,opts):
	bwa_args=[]
	if conf["bwa"]["aligner"]!="":
		bwa_args+=['mem']
	else:
		bwa_args+=[conf["bwa"]["aligner"]]
	if conf["bwa"]["threads"]!="":
		bwa_args+=['-t 8']
	else:
		bwa_args+=['-t '+ conf["bwa"]["threads"]]

	if conf["bwa"]["output-prefix"]!="":
		bwa_args+=[opts.out_path+'/BAM/' + conf["bwa"]["output-prefix"]]
	else:
		bwa_args+=[opts.out_path+'/BAM/' + opts.sample]

	path_bwa=conf["bwa"]["path"]

	return path_bwa,bwa_args



def Set_pirs_args(conf,opts):

	pirs_arguments=[]

	if conf["pirs"]["converage"]!="":
		pirs_arguments+=['--coverage='+conf["pirs"]["converage"]]
	else:
		pirs_arguments+=['--coverage=1000']
	if conf["pirs"]["read-len"]!="":
		pirs_arguments+=["--read-len="+conf["pirs"]["read-len"]]
	else:
		pirs_arguments+=["--read-len=150"]
	if conf["pirs"]["insert-len-mean"]!="":
		pirs_arguments+=["--insert-len-mean="+ conf["pirs"]["insert-len-mean"]]
	else:
		pirs_arguments+=["--insert-len-mean=300"]
	if conf["pirs"]["insert-len-sd"]!="":
		pirs_arguments+=["--insert-len-sd="+conf["pirs"]["insert-len-sd"]]
	else:
		pass
	if conf["pirs"]["base-calling-profile"]!="":
		pirs_arguments+=["--base-calling-profile="+conf["pirs"]["base-calling-profile"]]
	else:
		pirs_arguments+=["--no-subst-errors"]
	if conf["pirs"]["indel-error-profile"]!="":
		pirs_arguments+=["--indel-error-profile="+conf["pirs"]["indel-error-profile"]]
	else:
		pirs_arguments+=["--no-indel-errors"]
	if conf["pirs"]["gc-bias-profile"]!="":
		pirs_arguments+=["--gc-bias-profile="+conf["pirs"]["gc-bias-profile"]]
	else:
		pirs_arguments+=["--no-gc-content-bias"]
	if conf["pirs"]["error-rate"]!="":
		pirs_arguments+=["--error-rate="+conf["pirs"]["error-rate"]]
	else:
		pass
	if conf["pirs"]["substitution-error-algorithm"]!="":
		pirs_arguments+=["--substitution-error-algorithm="+conf["pirs"]["substitution-error-algorithm"]]
	else:
		pass
	if conf["pirs"]["quality-shift"]!="":
		pirs_arguments+=["--quality-shift="+conf["pirs"]["quality-shift"]]
	else:
		pass
	if conf["pirs"]["output-prefix"]!="":
		pirs_arguments+=["--output-prefix="+opts.out_path+'/FASTQ/' +conf["pirs"]["output-prefix"]]
	else:
		pirs_arguments+=["--output-prefix="+opts.out_path+'/FASTQ/' +opts.sample]
	if conf["pirs"]["threads"]!="":
		pirs_arguments+=["--threads="+conf["pirs"]["threads"]]
	else:
		pirs_arguments+=["--threads=4"]
	
	pirs_arguments+=conf["pirs"]["flags"].split(',')
	
	pirs_arguments+=["--output-file-type=gzip"]
	path_pirs=conf["pirs"]["path"]
	return path_pirs,pirs_arguments


def Set_picard(conf,opts):
	picard_args=[]
	if conf["picard"]["ram"]!="":
		picard_args+=[conf["picard"]["ram"]]
	else:
		picard_args+=[]

	path_picard=conf["picard"]["path"]
	return path_picard,picard_args


def Set_samtools(conf,opts):
	return conf["samtools"]["path"]

def Set_bedtools(conf,opts):
	return conf["bedtools"]["path"]

