# Frequently asked questions

Questions listed here are sorted by the executing order in the
[Running Guide](./RUNNING.md):

## About environments

#### What are system requirements?

The only system requirement is the size of the physical memory: 32G at least,
since the human genome indexes generated by STAR with default parameters will
require 30G memory to load when mapping reads.

#### Why to use container technology?

We introduced container technology here to assure the consistency of all running
environments and to speed up the booting of flexible cloud computing nodes.

#### Which version to use as reference?

At this very moment, GRCh38.p12 is the latest version of human genome reference.
However, current build of dbSNP (b151) is using GRCh38.p7 and current resource
bundle of GATK is based on GRCh38.p2.

Here, we choose to keep consistent with GATK:
- download genome reference, SNP references from GATK resource bundle;
- download genome annotations from [GENCODE](https://www.gencodegenes.org/releases/22.html).

## About alignment

#### Why to use STAR as the aligner?

We choose STAR due to [its robust performance with default parameters on
different situations](http://dx.doi.org/10.1038/nmeth.4106).

#### Why should we mark duplicates?

It is important in removing PCR duplicates which would introduce bias when
calling variants.

#### Who do we need base recalibration?

Base quality scores play an important role in weighing the evidence for or
against possible variant alleles during the variant discovery process, so it's
important to correct any systematic bias observed in the data.
