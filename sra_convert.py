#!/usr/bin/env python

from argparse import ArgumentParser
import os
import pypiper
import logmuse
import sys

__version__ = "0.1.0"

def _parse_cmdl(cmdl):
    parser = ArgumentParser(description="Automatic SRA data converter")
    # parser = pypiper.add_pypiper_args(parser, args=["output-parent"])
    parser.add_argument(
            "-m", "--mode", default="convert",
            choices=['convert', 'delete_sra', 'delete_bam', 'delete_fq'],
            help="What do you want to do? Default: convert")
    
    parser.add_argument(
            "-f", "--format", default="fastq", choices=['fastq', 'bam'],
            help="Convert to what format? Default: fastq")
    
    parser.add_argument(
            "-b", "--bamfolder", 
            default=safe_echo("SRABAM"),
            help="Optional: Specify a location to store bam files "
            "[Default: $SRABAM:" + safe_echo("SRABAM") + "]")

    parser.add_argument(
            "-q", "--fqfolder",
            default=safe_echo("SRAFQ"),
            help="Optional: Specify a location to store fastq files "
            "[Default: $SRAFQ:" + safe_echo("SRAFQ") + "]")
    
    parser.add_argument(
            "-s", "--srafolder", default=safe_echo("SRARAW"),
            help="Optional: Specify a location to store pipeline output "
            "[Default: $SRARAW:" + safe_echo("SRARAW") + "]")

    parser.add_argument(
            "--keep-sra", action='store_true', default=False,
            help="On convert mode, keep original sra data?")

    parser.add_argument( "-r", "--srr", required=True, nargs="+",
            help="SRR files")

    parser = pypiper.add_pypiper_args(parser, groups=["config", "logmuse"],
        args=["output-parent", "recover"])


    return parser.parse_args(cmdl)

def safe_echo(var):
    """ Returns an environment variable if it exists, or an empty string if not"""
    return os.getenv(var, "")

if __name__ == "__main__":
    args = _parse_cmdl(sys.argv[1:])

    global _LOGGER
    _LOGGER = logmuse.logger_via_cli(args)

    key=args.srr[0]

    if args.output_parent:
        outfolder = args.output_parent
    else:
        outfolder = os.path.join(args.srafolder, "sra_convert_pipeline")

    pm = pypiper.PipelineManager(name="sra_convert", outfolder=outfolder, args=args)

    nfiles = len(args.srr)
    failed_files = []

    for i in range(nfiles):
        srr_acc = os.path.splitext(os.path.basename(args.srr[i]))[0]
        pm.info("Processing {} of {} files: {}".format(str(i+1), str(nfiles), srr_acc))

        bamfile = os.path.join(args.bamfolder, srr_acc + ".bam")
        fq_prefix = os.path.join(args.fqfolder, srr_acc)
        
        if args.mode =="delete_sra" or not args.keep_sra:
            delete_sra = True
        else:
            delete_sra = False

        if args.mode == "convert":
            infile = args.srr[i]
            if not os.path.isfile(infile):
                pm.warning("Couldn't find sra file at: {}.".format(infile))
                failed_files.append(args.srr[i])
            if args.format == 'fastq':
                # fastq-dump --split-files will produce *_1.fastq and *_2.fastq
                # for paired-end data, and only *_1.fastq for single-end data.
                outfile = "{fq_prefix}_1.fastq.gz".format(fq_prefix=fq_prefix)
                cmd = "fastq-dump {data_source} --split-files --gzip -O {outfolder}".format(
                    data_source=infile, outfolder=args.fqfolder, nofail=True)
            elif args.format == 'bam':
                outfile = os.path.join(args.bamfolder, args.srr[i] + ".bam")
                cmd = "sam-dump -u {data_source} | samtools view -bS - > {outfile}".format(
                    data_source=infile, outfile=outfile, nofail=True)
            else:
                raise KeyError("Unknown format: {}".format(args.format))

            target = outfile
            ret = pm.run(cmd, target=target)
            if ret == 0:
                pm.info("Already completed files: {}" failed_files)
                try:
                    failed_files.remove(infile)
                except:
                    pass

        elif args.mode == "delete_bam":
            pm.timestamp("Cleaning bam file: {}".format(bamfile))
            pm.clean_add(bamfile)
        elif args.mode == "delete_fq":
            pm.timestamp("Cleaning fastq file(s): {}*".format(fq_prefix))
            fq_prefix = os.path.join(args.fqfolder, srr_acc)
            pm.clean_add("{fq_prefix}.fastq.gz".format(fq_prefix=fq_prefix))
            pm.clean_add("{fq_prefix}_[0-9].fastq.gz".format(fq_prefix=fq_prefix))
        if delete_sra:
            pm.timestamp("Cleaning sra file: {}".format(infile))
            pm.clean_add(infile)

    if len(failed_files) > 0:
        pm.fail_pipeline(Exception("Unable to locate the following files: {}".format(",".join(failed_files))))

    pm.stop_pipeline()
