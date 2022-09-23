import argparse
import polars as pl
import subprocess
import os
import sys

parser = argparse.ArgumentParser(prog="regenie_lamp.py",
                                 description="Python wrapper for regenie GWAS analysis using plink files. \
                                 Performs steps 1 and 2, as well as generates QQ and manhattan plots.")
i_o_group = parser.add_argument_group("Input/Output")
i_o_group.add_argument("--regenie-binary", type=str, help="path to regenie binary, if in PATH put path name of command")
i_o_group.add_argument("--out-dir", type=str, help="path to output directory")
i_o_group.add_argument("--out-prefix", type=str, help="prefix for output files")
i_o_group.add_argument("--threads", type=int, help="number of threads to use (default: 8)", default=8)
# Genetic Data Set Arguments
genetics_group = parser.add_argument_group("Genetic Data Files - Only specify either --pfile or --bfile")
genetics_group.add_argument("--pfile", help="Plink2 PGEN file prefix. All PGEN files must have the same prefix", default=None)
genetics_group.add_argument("--bfile", help="Plink BED file prefix. All BED files must have the same prefix", default=None)
genetics_group.add_argument("--step-one-variants", help="File containing variants to be used in step 1.", default=None)
genetics_group.add_argument("--step-two-variants", help="File containing variants to be used in step 2.", default=None)
# Sample Data Arguments
sample_data = parser.add_argument_group("Sample Data Files - Assumes Tab Separated")
sample_data.add_argument("--phenotypes",
                         help="Phenotype file, must contain either all binary or all continuous phenotypes")
sample_data.add_argument("--bt", help="Flag to indicate that the phenotype file contains binary traits", action="store_true")
sample_data.add_argument("--covariates", help="Covariate file, rows with missing records will be dropped", default=None)
sample_data.add_argument("--sample-ids", help="File containing sample IDs to be used in analysis", default=None)

args = parser.parse_args()


def determine_genetic_files(args):
    if args.pfile is None and args.bfile is None:
        sys.exit("Error: No genetic data set provided. Please provide either a PGEN or BED file prefix.")
    if args.pfile is not None and args.bfile is not None:
        raise ValueError("Please specify either --pfile or --bfile, not both")
    elif args.pfile is not None:
        prefix = args.pfile
        genetics = prefix + ".pgen"
        variants = prefix + ".pvar"
        samples = prefix + ".psam"
    elif args.bfile is not None:
        prefix = args.bfile
        genetics = prefix + ".bed"
        variants = prefix + ".bim"
        samples = prefix + ".fam"

    for file in [genetics, variants, samples]:
        if not os.path.isfile(file):
            sys.exit("Error: {} does not exist".format(file))
    return prefix, genetics, variants, samples


def compare_sample_ids(genetic_samples, args):
    samples = pl.read_csv(genetic_samples, sep="\t")
    if args.phenotypes is not None:
        if os.path.isfile(args.phenotypes):
            phenotypes = pl.read_csv(args.phenotypes, sep="\t", null_values='NA').rename({"FID": "#FID"})
            if samples.join(phenotypes, on=["#FID", "IID"]).shape[0] == 0:
                sys.exit("Error: No samples in common between genetic data and phenotypes file")
        else:
            sys.exit("Error: Phenotype file {} does not exist".format(args.phenotypes))
    else:
        sys.exit("Error: A phenotype file must be provided")

    if args.sample_ids is not None:
        if os.path.isfile(args.sample_ids):
            sample_ids = pl.read_csv(args.sample_ids, sep="\t", new_columns=["#FID", "IID"])
            if len(sample_ids.columns) > 2:
                sys.exit("Error: sample_ids file must contain only two columns")
            if samples.join(sample_ids, on=["#FID", "IID"]).shape[0] == 0:
                sys.exit("Error: No samples in common between genetic data and sample_ids file")
        else:
            sys.exit("Error: Sample ID file {} does not exist".format(args.sample_ids))

    if args.covariates is not None:
        if os.path.isfile(args.covariates):
            covariates = pl.read_csv(args.covariates, sep="\t", null_values='NA').rename({"FID": "#FID"})
            covariates = covariates.drop_nulls()
            if samples.join(covariates, on=["#FID", "IID"]).shape[0] == 0:
                sys.exit("Error: No samples in common between genetic data and covariates file")
        else:
            sys.exit("Error: Covariates file {} does not exist".format(args.covariates))


def compare_variant_ids(genetic_variants, args):
    variants = pl.read_csv(genetic_variants, sep="\t", dtypes={"#CHROM": pl.Utf8})
    if args.step_one_variants is not None:
        if os.path.isfile(args.step_one_variants):
            step_one_variants = pl.read_csv(args.step_one_variants, sep="\t", new_columns=['ID'])
            if variants.join(step_one_variants, on=["ID"]).shape[0] == 0:
                sys.exit("Error: No variants in common between genetic data and step_one_variants file")
        else:
            sys.exit("Error: Step One Variants file {} does not exist".format(args.step_one_variants))
    if args.step_two_variants is not None:
        if os.path.isfile(args.step_two_variants):
            step_two_variants = pl.read_csv(args.step_two_variants, sep="\t", new_columns=['ID'])
            if variants.join(step_two_variants, on=["ID"]).shape[0] == 0:
                sys.exit("Error: No variants in common between genetic data and step_two_variants file")
        else:
            sys.exit("Error: Step Two Variants file {} does not exist".format(args.step_two_variants))


def main(args):
    prefix, genetics, variants, samples = determine_genetic_files(args)
    compare_sample_ids(samples, args)
    compare_variant_ids(variants, args)
    for folder in ['step_one_out', 'step_two_out', f'tmp', 'plots']:
        if not os.path.exists(f'{args.out_dir}/{folder}'):
            os.makedirs(f'{args.out_dir}/{folder}')

    step_one_command = f"""{args.regenie_binary} \
    --step 1 \
    --bsize 1000 \
    --threads {args.threads}
    --lowmem \
    --lowmem-prefix {args.out_dir}/tmp/tmp_{prefix} \
    {f"--bed {prefix}" if args.bfile is not None else f"--pgen {prefix}"} \
    --extract {args.step_one_variants} \
    --phenoFile {args.phenotypes} {"--bt" if args.bt else ""} \
    {f"--covarFile {args.covariates}" if args.covariates is not None else ""} \
    {f"--keep {args.sample_ids}" if args.sample_ids is not None else ""} \
    --out {args.out_dir}/step_one_out/{args.out_prefix}"""

    subprocess.run(step_one_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    step_two_command = f"""{args.regenie_binary} \
    --step 2 \
    --bsize 1000 \
    --threads {args.threads} \
    {f"--bed {prefix}" if args.bfile is not None else f"--pgen {prefix}"} \
    {f"--extract {args.step_two_variants}" if args.step_two_variants is not None else ""} \
    --phenoFile {args.phenotypes} {"--bt" if args.bt else ""} \
    {f"--covarFile {args.covariates}" if args.covariates is not None else ""} \
    {f"--keep {args.sample_ids}" if args.sample_ids is not None else ""} \
    --pred {args.out_dir}/step_one_out/{args.out_prefix}_pred.list \
    --firth --approx \
    --out {args.out_dir}/step_two_out/{args.out_prefix}"""

    subprocess.run(step_two_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


if __name__ == "__main__":
    main(args)
