#!/usr/bin/env python

"""
Simulate crossing each sample in a VCF file against the reference sequence.

For each call in the file:
- homozygous ref sites stay homozygous ref
- homozygous alt sites become heterozygous
- all other sites become unknown (no call)
"""

import sys
import vcf
import numpy
import click
import collections

# Location of ss deficiency in dm6 coordinates
def_chrom = "3R"
def_start = 16305713
def_end = 16502732

class Stats( object ):
    def __init__( self ):
        self.non_snp = self.multi_allelic = self.missing = self.ref = self.passed = 0

class B( bytes ):
    def encode( self ):
        return self

@click.command()
@click.argument( 'infname' )
@click.argument( 'outfname' )
def main( infname, outfname ):
    # Parse VCF from stdin -- htslib requires a filename
    vcf_reader = vcf.VCFReader( open( infname ) )
    out = open( outfname, "w" )
    # Start new VCF on stdout
    out.write( "##fileformat=VCFv4.1\n" )
    out.write( "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join( vcf_reader.samples ) + "\n" )
    # Track some stats
    stats = Stats()
    # Visit each variant
    for variant in vcf_reader:
        # Ignore indels or SVs
        if not ( variant.is_snp or variant.is_indel ):
            stats.non_snp += 1
            continue
        # Ignore anything that is not bi-allelic
        if len( variant.ALT ) > 1:
            stats.multi_allelic += 1
            continue
        # "Cross" each individual against reference
        #   One-based inclusive per GFF and VCF formats
        gt_types = [ 3 if call.gt_type is None else call.gt_type for call in variant ]
        if variant.CHROM == def_chrom and def_start <= variant.POS < def_end:
            new_gt_types = cross_gt_types( gt_types, homozygous_becomes=2 )
        else:
            new_gt_types = cross_gt_types( gt_types )
        # If all missing or reference after cross, skip site
        if ( new_gt_types == 3 ).all():
            stats.missing += 1
            continue
        if ( ( new_gt_types == 0 ) | ( new_gt_types == 3 ) ).all():
            stats.ref += 1
            continue
        write_variant( out, variant, new_gt_types )
        stats.passed += 1
    # Done
    sys.stderr.write( "Dropped:" )
    sys.stderr.write( "  %d non snp" % stats.non_snp )
    sys.stderr.write( "  %d multi allelic" % stats.multi_allelic )
    sys.stderr.write( "  %d missing" % stats.missing )
    sys.stderr.write( "  %d all reference" % stats.ref )
    sys.stderr.write( "Wrote:" )
    sys.stderr.write( "  %d variants" % stats.passed )

def cross_gt_types( gt_types, homozygous_becomes=1 ):
    """
    Given an array of genotype types (0, 1, 2, 3: missing), simulate crossing
    each individual against the reference. Site that are heterozygous will
    become missing.
    """
    new_gt_types = numpy.array( gt_types )
    n_missing = 0
    for i in range( len( new_gt_types ) ):
        gt_type = new_gt_types[i]
        if gt_type == 0 or gt_type == 3:
            # Missing calls and homozygous reference calls are unchanged
            continue
        elif gt_type == 2:
            # Homozygous alt becomes heterozygous
            gt_type = 1
        else:
            # Heterozygous sites become missing
            gt_type = 3
        # Update array with new type
        new_gt_types[i] = gt_type
    return new_gt_types

def write_variant( out, variant, new_gt_types ):
    """
    Write a minimal variant record for `variant` using genotypes from
    `new_gt_types`.
    """
    out.write( "\t".join( map( str, (
        variant.CHROM, variant.POS, variant.ID,
        variant.REF, variant.ALT[0], ".", ".", ".", "GT" ) ) ) )
    for gt_type in new_gt_types:
        out.write( "\t" )
        if gt_type == 0:
            out.write( "0" )
        elif gt_type == 1:
            out.write( "1" )
        elif gt_type == 1:
            out.write( "2" )
        else:
            out.write( "." )
    out.write( "\n" )

if __name__ == "__main__":
    main()
